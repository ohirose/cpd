// Copyright (c) 2015 Osamu Hirose
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

/*-- NOTE --------------------------------------------------------------------o
| Simple implementation of Coherent Point Drift, a method of point set        |
| registration. This implementation uses "dposv" in Lapack for solving linear |
| equations. Be careful when linking lapack library (dependent on OS).        |
o-----------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<string.h>
#include<math.h>
#include"util.h"
#include"func.h"

double ** readPoints(int *nr, int *nc, const char *file, const char mode){
  int n,d,L,N,D,si,sd; FILE *fp; double **X,*b;
  si=sizeof(int   );
  sd=sizeof(double);

  fp=fopen(file,"rb");if(!fp){printf("File not found: %s\n",file);exit(EXIT_FAILURE);}
  fread(nr,si,1,fp);
  fread(nc,si,1,fp); N=*nr;D=*nc;L=N*D;

  if(N<D){
    printf("The number of points in a point set is at least \n");
    printf("the dimension of the vector space. Abort.       \n"); exit(1);
  }

  X=calloc2d(N,D);b=malloc(sd*L);fread(b,sd,L,fp);fclose(fp);
  switch(mode){
    case 'r': for(n=0;n<N;n++)for(d=0;d<D;d++) X[n][d]=b[d+n*D]; break;
    case 'c': for(n=0;n<N;n++)for(d=0;d<D;d++) X[n][d]=b[n+d*N]; break;
  } free(b);

  return X;
}

int writePoints(const char *file, const double **X, const int N, const int D){
  int n,d; FILE *fp=fopen(file,"w");if(!fp){printf("Can't open: %s\n",file);exit(EXIT_FAILURE);}
  for(n=0;n<N;n++)for(d=0;d<D;d++)
    fprintf(fp,"%lf%c",X[n][d],d==D-1?'\n':'\t');
  fclose(fp);
  return 0;
}

int scalePoints(double **X, const int N, const int D, const double *scale){
  int n,d; for(n=0;n<N;n++)for(d=0;d<D;d++) X[n][d]*=scale[d];
  return 0;
}

int normPoints(double **X, double *mu, double *sgm, const int N, const int D){
  int n,d; double val=0;
  for(d=0;d<D;d++){mu[d]=0;for(n=0;n<N;n++) mu[d]+=X[n][d];mu[d]/=N;}
  for(n=0;n<N;n++)for(d=0;d<D;d++) val +=SQ(X[n][d]-mu[d]); val/=N*D;*sgm=sqrt(val);
  for(n=0;n<N;n++)for(d=0;d<D;d++) X[n][d]=(X[n][d]-mu[d])/(*sgm);
  return 0;
}

int revertPoints(double **X, double *mu, double *sgm, const int N, const int D){
  int n,d; for(n=0;n<N;n++)for(d=0;d<D;d++) X[n][d]=X[n][d]*(*sgm)+mu[d];
  return 0;
}

int readPrms(double prms[5],const char *file){
  FILE* fp=fopen(file,"r");if(!fp){printf("File not found: %s\n",file);exit(EXIT_FAILURE);}
  fscanf(fp,"nlp:%lf\n", prms+0);
  fscanf(fp,"omg:%lf\n", prms+1);
  fscanf(fp,"bet:%lf\n", prms+2);
  fscanf(fp,"lmd:%lf\n", prms+3);
  fscanf(fp,"dz:%lf\n",  prms+4);
  fclose(fp);
  return 0;
}

int printOptProcess(const char *file, const double *Q, const int lp, const int M, const int D){
  if(Q){FILE *fp=fopen(file,"wb");if(!fp){printf("Can't open: %s\n",file);exit(EXIT_FAILURE);}
    fwrite(&M, sizeof(int),   1,     fp);
    fwrite(&D, sizeof(int),   1,     fp);
    fwrite(&lp,sizeof(int),   1,     fp);
    fwrite(Q,  sizeof(double),lp*M*D,fp);
    fclose(fp);
  } return 0;
}

int main(int argc, char **argv){
  int M,N,D,nlp,nlpr[3]={0,0,0},size[3],flag=0,verb; double prms[6];
  double **W,**T,**G,**P,***C,*A,*B,**X,**Y,**Z,*Q0,*Q1,*Q2;
  double sgmX,sgmY,muX[3],muY[3],asp[3]={1,1,1},iasp[3]={1,1,1};

  if(argc!=4){
    printf("./cpd <mode> <X> <Y>\n\n");
    printf("NOTE: At least one of characters r, a, and c must be included in <mode>.       \n");
    printf("Optionally, m and i can be attatched to <mode> which specifies a print option. \n");
    printf("r: rotate, a: affine, c: cpd, i: information, m: memorize optimization process.\n");
    printf("\n");
    exit(1);
  }

  flag+=strchr(argv[1],(int)'r')!=NULL?1:0;
  flag+=strchr(argv[1],(int)'a')!=NULL?2:0;
  flag+=strchr(argv[1],(int)'c')!=NULL?4:0;
  flag+=strchr(argv[1],(int)'m')!=NULL?8:0;
  verb =strchr(argv[1],(int)'i')!=NULL?1:0;

  readPrms(prms,"prms.txt");
  X=readPoints(&N,&D,argv[2],'r');
  Y=readPoints(&M,&D,argv[3],'r');
  size[0]=M;size[1]=N;size[2]=D;nlp=prms[0];

  W = calloc2d(M,  D  ); A = calloc  (M*M,sizeof(double));
  T = calloc2d(M,  D  ); B = calloc  (M*D,sizeof(double));
  G = calloc2d(M,  M  ); C = calloc3d(4,N>M?N:M,D);
  P = calloc2d(M+1,N+1);

  Q0=flag&1?malloc(nlp*M*D*sizeof(double)):NULL;
  Q1=flag&2?malloc(nlp*M*D*sizeof(double)):NULL;
  Q2=flag&4?malloc(nlp*M*D*sizeof(double)):NULL;

  #define CD  (const double **)
  #define CD1 (const double * )
  if(D==3) {asp[2]=prms[4];iasp[2]=1.0/prms[4];}
  scalePoints(X,N,D,asp); normPoints(X,muX,&sgmX,N,D);
  scalePoints(Y,M,D,asp); normPoints(Y,muY,&sgmY,M,D);

  if(flag&1) {nlpr[0]=rot   (W,T,P,C,Q0,CD X,CD Y,size,prms,verb); if(flag&6){Z=Y;Y=T;T=Z;}}
  if(flag&2) {nlpr[1]=affine(W,T,P,C,Q1,CD X,CD Y,size,prms,verb); if(flag&4){Z=Y;Y=T;T=Z;}}
  if(flag&4) {nlpr[2]=cpd   (W,T,G,P,C[0],A,B,Q2,CD X,CD Y,size,prms,verb);}

  writePoints("T1.txt", CD T,M,D);
  writePoints("X1.txt", CD X,N,D);
  writePoints("Y1.txt", CD Y,M,D);

  if(Q0) printOptProcess("otw-r.bin",CD1 Q0,nlpr[0],M,D);
  if(Q1) printOptProcess("otw-a.bin",CD1 Q1,nlpr[1],M,D);
  if(Q2) printOptProcess("otw-c.bin",CD1 Q2,nlpr[2],M,D);

  return 0;
}

