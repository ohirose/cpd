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
| registration. This implementation uses "dposv", "dsyev", and "dgesdd" in    |
| Lapack. Be careful when linking lapack library (dependent on OS).           |
o-----------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<string.h>
#include<math.h>
#include<unistd.h>
#include"util.h"
#include"func.h"

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

int readPrms(double prms[6],const char *file){
  FILE* fp=fopen(file,"r");if(!fp){printf("File not found: %s\n",file);exit(EXIT_FAILURE);}
  fscanf(fp,"nlp:%lf\n", prms+0);
  fscanf(fp,"omg:%lf\n", prms+1);
  fscanf(fp,"bet:%lf\n", prms+2);
  fscanf(fp,"lmd:%lf\n", prms+3);
  fscanf(fp,"rnk:%lf\n", prms+4);
  fscanf(fp,"dz:%lf\n",  prms+5);
  fclose(fp);
  return 0;
}

int checkPrms(double prms[6]){ int nlp,K; double omg,lmd,bet,dz;
  nlp=prms[0];omg=prms[1];lmd=prms[2];bet=prms[3];K=prms[4];dz=prms[5];
  if(!(nlp>0))       {printf("-n: Argument must be a positive integer. Abort.\n");        exit(EXIT_FAILURE);}
  if(!(omg>0&&omg<1)){printf("-w: Argument must range (0,1). Abort.\n");                  exit(EXIT_FAILURE);}
  if(!(lmd>0))       {printf("-l: Argument must be positive. Abort.\n");                  exit(EXIT_FAILURE);}
  if(!(bet>0))       {printf("-b: Argument must be positive. Abort.\n");                  exit(EXIT_FAILURE);}
  if(!(dz >0))       {printf("-z: Argument must be positive. Abort.\n");                  exit(EXIT_FAILURE);}
  if(!(K >=0))       {printf("-r: Argument must be a positive integer or zero. Abort.\n");exit(EXIT_FAILURE);}
  return 0;
}

int printOptProcess(const char *file, const double *S, const int lp, const int M, const int D){
  if(S){FILE *fp=fopen(file,"wb");if(!fp){printf("Can't open: %s\n",file);exit(EXIT_FAILURE);}
    fwrite(&M, sizeof(int),   1,     fp);
    fwrite(&D, sizeof(int),   1,     fp);
    fwrite(&lp,sizeof(int),   1,     fp);
    fwrite(S,  sizeof(double),lp*M*D,fp);
    fclose(fp);
  } return 0;
}

int main(int argc, char **argv){
  int K,M,N,D,nlp,nlpr[3]={0,0,0},size[3],flag=0,verb,optc;
  double **W,**T,**G,**P,***C,***U,***V,*A,*B,**X,**Y,**Z,*S0,*S1,*S2;
  double dz,sgmX,sgmY,muX[3],muY[3],asp[3]={1,1,1},iasp[3]={1,1,1};
  char mode,opt,**optp,*fout,ftxt[32]="T.txt",fbin[32]="T.bin",fprm[32]="";
  double prms[6]={1000,0.1,1.0,1.0,0,1.0};/*default: nloop,omega,lambda,beta,rank,zscale*/

  if(argc<4){
    printf("\n");
    printf("  USAGE: ./cpd <mode> <X> <Y> (+options) \n\n");
    printf("  MODE: At least one of characters r, a, and c must be included in <mode>.           \n");
    printf("  Optionally, m and i which specify print options can be attatched to <mode>.        \n");
    printf("  r: rotate, a: affine, c: cpd, i: information, m: memorize optimization process.  \n\n");
    printf("  OPTIONs: options must be added AFTER the arguments.                                \n");
    printf("  -n: nloop, -w omega, -l lambda, -b beta, -r rank, -z zscale, -p <parameter file>.\n\n");
    printf("  EXAMPLE: ./cpd raci X.txt Y.txt -w 0.5 -l 15 -b 0.9 -z 3.5 -n 2000               \n\n");
    exit(1);
  }

  flag+=strchr(argv[1],(int)'r')!=NULL?1:0;
  flag+=strchr(argv[1],(int)'a')!=NULL?2:0;
  flag+=strchr(argv[1],(int)'c')!=NULL?4:0;
  flag+=strchr(argv[1],(int)'m')!=NULL?8:0;
  verb =strchr(argv[1],(int)'i')!=NULL?1:0;

  optp=argv+3; optc=argc-3;
  while((opt=getopt(optc,optp,"n:w:l:b:r:z:p:"))!=-1){
    switch(opt){
      case 'n': prms[0]=atof(optarg);  break;
      case 'w': prms[1]=atof(optarg);  break;
      case 'l': prms[2]=atof(optarg);  break;
      case 'b': prms[3]=atof(optarg);  break;
      case 'r': prms[4]=atof(optarg);  break;
      case 'z': prms[5]=atof(optarg);  break;
      case 'p': strcpy (fprm,optarg);  break;
      default : exit(EXIT_FAILURE);
    }
  } if(strlen(fprm)) readPrms(prms,fprm); checkPrms(prms);

  X=read2d(&N,&D,&mode,argv[2]);
  Y=read2d(&M,&D,&mode,argv[3]);
  size[0]=M;size[1]=N;size[2]=D;nlp=prms[0];K=prms[4];dz=prms[5];
  fout=(mode=='t')?ftxt:fbin;

  if(verb){
    printf("--------------- SUMMARY ---------------------------------------------\n");
    printf("\n");
    printf("  Input Data:\n");
    printf("    Point set 1 (reference): [%s]\n", argv[2]);
    printf("    Point set 2 (floating):  [%s]\n", argv[3]);
    printf("    Size of point set 1: [%3d,%2d]\n", M,D);
    printf("    Size of point set 2: [%3d,%2d]\n", N,D); printf("\n");
    printf("  Deformation Model:\n    [");
    if(flag&1) printf("rigid" ); if(flag&1&&(flag&2||flag&4)) printf("+");
    if(flag&2) printf("affine"); if(flag&2&&(flag&4))         printf("+");
    if(flag&4) printf("cpd"   ); printf("]\n\n");
    printf("  Parameters: ");if(strlen(fprm))printf("%s",fprm); printf("\n");
    printf("    nloop  = %d\n", (int) prms[0]);
    printf("    omega  = %lf\n",      prms[1]);
    printf("    lambda = %lf\n",      prms[2]);
    printf("    beta   = %lf\n",      prms[3]);
    if(K)    printf("    rank   = %d\n", (int) prms[4]);
    if(D==3) printf("    zscale = %lf\n",prms[5]);
    if(strlen(fprm)){printf("\n");
      printf("  ** Parameters specified by %s were used, that is,\n", fprm);
      printf("  those specified by command-line arguments are disabled.\n");
    }
    printf("\n");
    printf("  Output: \n");
    printf("    [%s]\n\n",fout);
  }

  W = calloc2d(M,  D  ); A = calloc  (M*M,sizeof(double));
  T = calloc2d(M,  D  ); B = calloc  (M*M,sizeof(double));
  G = calloc2d(M,  M  ); C = calloc3d(4,N>M?N:M,D);
  U = calloc3d(2,K+1,M); V = calloc3d(2,M,D);
  P = calloc2d(M+1,N+1);

  S0=flag&1?malloc(nlp*M*D*sizeof(double)):NULL;
  S1=flag&2?malloc(nlp*M*D*sizeof(double)):NULL;
  S2=flag&4?malloc(nlp*M*D*sizeof(double)):NULL;

  #define CD  (const double **)
  #define CD1 (const double * )
  if(D==3) {asp[2]=dz;iasp[2]=1.0/dz;}
  scalePoints(X,N,D,asp); normPoints(X,muX,&sgmX,N,D);
  scalePoints(Y,M,D,asp); normPoints(Y,muY,&sgmY,M,D);

  if(flag&1) {nlpr[0]=rot   (W,T,P,C,S0,CD X,CD Y,size,prms,verb); if(flag&6){Z=Y;Y=T;T=Z;}}
  if(flag&2) {nlpr[1]=affine(W,T,P,C,S1,CD X,CD Y,size,prms,verb); if(flag&4){Z=Y;Y=T;T=Z;}}
  if(flag&4) {nlpr[2]=cpd   (W,T,P,G,U,V,A,B,S2,CD X,CD Y,size,prms,verb);}

  revertPoints(T,muX,&sgmX,M,D); scalePoints(T,M,D,iasp);
  revertPoints(Y,muY,&sgmY,M,D); scalePoints(Y,M,D,iasp);
  revertPoints(X,muX,&sgmX,N,D); scalePoints(X,N,D,iasp);

  write2d(fout,CD T,M,D);

  if(S0) printOptProcess("otw-r.bin",CD1 S0,nlpr[0],M,D);
  if(S1) printOptProcess("otw-a.bin",CD1 S1,nlpr[1],M,D);
  if(S2) printOptProcess("otw-c.bin",CD1 S2,nlpr[2],M,D);

  if(verb) printf("---------------------------------------------------------------------\n");

  return 0;
}

