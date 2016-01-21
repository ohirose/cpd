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
#include<ctype.h>
#include"util.h"
#include"func.h"
#include"info.h"

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

void readPrms(double prms[6],const char *file){
  int i,l,L,N=6; FILE *fp; char *s, line[1024],*names[]={"nloop","omega","lambda","beta","rank","zscale"};

  fp=fopen(file,"r");if(!fp){printf("File not found: %s\n",file);exit(EXIT_FAILURE);}
  while(fgets(line,1024,fp)){
    if(line[0]=='#') continue; L=strlen(line);
    for(l=0;l<L;l++) line[l]=tolower(line[l]);
    for(i=0;i<N;i++){s=names[i];if(strstr(line,s))sscanf(line+strlen(s),"%lf\n",&prms[i]);}
  } fclose(fp);

  return;
}

int checkPrms(double prms[6]){ int nlp,rank; double omg,lmd,bet,dz;
  nlp=prms[0];omg=prms[1];lmd=prms[2];bet=prms[3];rank=prms[4];dz=prms[5];
  if(!(nlp>0))       {printf("-n: Error: Argument must be a positive integer. Abort.\n"); exit(EXIT_FAILURE);}
  if(!(omg>0&&omg<1)){printf("-w: Error: Argument must be in range (0,1). Abort.\n");     exit(EXIT_FAILURE);}
  if(!(lmd>0))       {printf("-l: Error: Argument must be positive. Abort.\n");           exit(EXIT_FAILURE);}
  if(!(bet>0))       {printf("-b: Error: Argument must be positive. Abort.\n");           exit(EXIT_FAILURE);}
  if(!(dz >0))       {printf("-z: Error: Argument must be positive. Abort.\n");           exit(EXIT_FAILURE);}
  if(!(rank>=0))     {printf("-r: Error: Argument must be a positive integer. Abort.\n"); exit(EXIT_FAILURE);}
  return 0;
}

int main(int argc, char **argv){
  int rank,M,N,D,nlp,nlpr[3]={0,0,0},size[3],flag=0,verb,optc,sd=sizeof(double);
  double **W,**T,**G,**P,***C,***U,***V,*A,*B,**X,**Y,**Z,*S0,*S1,*S2;
  double dz,sgmX,sgmY,muX[3],muY[3],asp[3]={1,1,1},iasp[3]={1,1,1};
  char mode,opt,**optp,*files[4],fout[256]="",fprm[256]="",ftxt[32]="T.txt",fbin[32]="T.bin";
  double prms[6]={1000,0.1,1.0,1.0,0,1.0};/*default: nloop,omega,lambda,beta,rank,zscale*/

  if      (argc==2){printVersion(); exit(EXIT_SUCCESS);}
  else if (argc <4){printUsage  (); exit(EXIT_SUCCESS);}

  flag+=strchr(argv[1],(int)'r')!=NULL?1:0;
  flag+=strchr(argv[1],(int)'a')!=NULL?2:0;
  flag+=strchr(argv[1],(int)'c')!=NULL?4:0;
  flag+=strchr(argv[1],(int)'m')!=NULL?8:0;
  verb =strchr(argv[1],(int)'q')!=NULL?0:1;

  optp=argv+3; optc=argc-3;
  while((opt=getopt(optc,optp,"n:w:l:b:r:z:p:o:v"))!=-1){
    switch(opt){
      case 'n': prms[0]=atof(optarg);  break;
      case 'w': prms[1]=atof(optarg);  break;
      case 'l': prms[2]=atof(optarg);  break;
      case 'b': prms[3]=atof(optarg);  break;
      case 'r': prms[4]=atof(optarg);  break;
      case 'z': prms[5]=atof(optarg);  break;
      case 'p': strcpy (fprm,optarg);  break;
      case 'o': strcpy (fout,optarg);  break;
      case 'v': printVersion();        break;
      default : exit(EXIT_FAILURE);
    }
  }

  X=read2d(&N,&D,&mode,argv[2]);
  Y=read2d(&M,&D,&mode,argv[3]);

  if(D< 2&&D> 3){printf("Error: Dimension must be 2 or 3.             \n");exit(EXIT_FAILURE);}
  if(N<=D||M<=D){printf("Error: #points must be greater than dimension\n");exit(EXIT_FAILURE);}

  if( strlen(fprm)) readPrms(prms,fprm); checkPrms(prms);
  if(!strlen(fout)) strcpy(fout,(mode=='t')?ftxt:fbin);

  files[0]=argv[2]; size[0]=M; nlp =prms[0];
  files[1]=argv[3]; size[1]=N; rank=prms[4];
  files[2]=fprm;    size[2]=D; dz  =prms[5];
  files[3]=fout;

  if(verb) printInfo((const char**)files,size,prms,flag);

  W = calloc2d(M,D);        A = calloc  (M*M,sd);
  T = calloc2d(M,D);        B = calloc  (M*M,sd);
  G = calloc2d(M,M);        C = calloc3d(4,N>M?N:M,D);
  U = calloc3d(2,rank+1,M); P = calloc2d(M+1,N+1);
  V = calloc3d(2,M,D);

  S0=(flag&1 && flag&8)?malloc(nlp*M*D*sd):NULL;
  S1=(flag&2 && flag&8)?malloc(nlp*M*D*sd):NULL;
  S2=(flag&4 && flag&8)?malloc(nlp*M*D*sd):NULL;

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

  return 0;
}

