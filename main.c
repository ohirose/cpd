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

  fp=fopen(file,"rb");
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
  int n,d; FILE *fp=fopen(file,"w");
  for(n=0;n<N;n++)for(d=0;d<D;d++)
    fprintf(fp,"%lf%c",X[n][d],d==D-1?'\n':'\t');
  fclose(fp);
  return 0;
}

int normPoints(double **X, const int N, const int D){
  int n,d; double val,sgm;
  for(d=0;d<D;d++){val=sgm=0;
    for(n=0;n<N;n++)val+=X[n][d];val/=N;
    for(n=0;n<N;n++)X[n][d]-=val;
  }
  for(n=0;n<N;n++)for(d=0;d<D;d++)sgm+=SQ(X[n][d]);sgm/=N*D;sgm=sqrt(sgm);
  for(n=0;n<N;n++)for(d=0;d<D;d++)X[n][d]/=sgm;
  return 0;
}

int rescalePoints(double **X, const int N, const int D, const double dz){
  int i; if(D==3&&fabs(dz-1.0)>1e-8) for(i=0;i<N;i++) X[i][2]*=dz;
  return 0;
}

int readPrms(double prms[5],const char *file){
  FILE* fp=fopen(file,"r");
  fscanf(fp,"nlp:%lf\n", prms+0);
  fscanf(fp,"omg:%lf\n", prms+1);
  fscanf(fp,"bet:%lf\n", prms+2);
  fscanf(fp,"lmd:%lf\n", prms+3);
  fscanf(fp,"dz:%lf\n",  prms+4);
  fclose(fp);
  return 0;
}

int printOptProcess(const char *file, const double *Q, const int lp, const int M, const int D){
  if(Q){FILE *fp=fopen(file,"wb");
    fwrite(&M, sizeof(int),   1,     fp);
    fwrite(&D, sizeof(int),   1,     fp);
    fwrite(&lp,sizeof(int),   1,     fp);
    fwrite(Q,  sizeof(double),lp*M*D,fp);
    fclose(fp);
  } return 0;
}

int main(int argc, char **argv){
  int M,N,D,nlp,nlpr[3],size[3],flag=0,verb;double prms[6];
  double **W,**T,**G,**P,***C,*A,*B,**X,**Y,**Z,*Q0,*Q1,*Q2;

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

  rescalePoints (X,N,D,prms[4]);
  rescalePoints (Y,M,D,prms[4]);
  normPoints    (X,N,D);
  normPoints    (Y,M,D);

  #define CD (const double **)
  if(flag&1) {nlpr[0]=rot   (W,T,P,C,Q0,CD X,CD Y,size,prms,verb); if(flag&6){Z=Y;Y=T;T=Z;}} 
  if(flag&2) {nlpr[1]=affine(W,T,P,C,Q1,CD X,CD Y,size,prms,verb); if(flag&4){Z=Y;Y=T;T=Z;}}
  if(flag&4) {nlpr[2]=cpd   (W,T,G,P,C[0],A,B,Q2,CD X,CD Y,size,prms,verb);}

  writePoints("T1.txt", CD T,M,D);
  writePoints("X1.txt", CD X,N,D);
  writePoints("Y1.txt", CD Y,M,D);
  #undef  CD

  #define CD (const double * )
  if(Q0) printOptProcess("otw-r.bin",CD Q0,nlpr[0],M,D);
  if(Q1) printOptProcess("otw-a.bin",CD Q1,nlpr[1],M,D);
  if(Q2) printOptProcess("otw-c.bin",CD Q2,nlpr[2],M,D);
  #undef  CD1

  return 0;
}

