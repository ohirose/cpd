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
#include<math.h>
#include"util.h"

int cpd(double       **  W,        /*  M   x  D   | displacement matrix   */
        double       **  T,        /*  M   x  D   | Moved reference       */
        double       **  G,        /*  M   x  M   | Gram matrix of Y      */
        double       **  P,        /*  M+1 x  N+1 | Assighment probablity */
        double       **  C,        /*  M   x  D   | Working wemory (2D)   */
        double       *   A,        /*  M   x  M   | Working wemory (1D)   */
        double       *   B,        /*  M   x  D   | Working wemory (1D)   */
        const double **  X,        /*  N   x  D   | Point set 1 (Data)    */
        const double **  Y,        /*  N   x  D   | Point set 2 (Data)    */
        const int        size[3],  /*  M,  N, D                           */
        const double     prms[4]   /*  omg,bet,lmd,nloop                  */
       ){
  int    i,m,n,d,M,N,D,lp,nloop,info; char uplo='U';
  double omg,lmd,bet2,reg=1e-9;
  double sgm2=0,noise,val,pres1=1e10,pres2=1e20,pren,conv;

  M=size[0];N=size[1];D=size[2];
  omg=prms[0];bet2=SQ(prms[1]);lmd=prms[2];nloop=(int)prms[3];

  /* initialize */
  for(m=0;m<M;m++)for(d=0;d<D;d++) W[m][d]=0;
  for(m=0;m<M;m++)for(d=0;d<D;d++) T[m][d]=Y[m][d];
  for(m=0;m<M;m++)for(n=0;n<N;n++) sgm2+=dist2(X[n],Y[m],D);sgm2/=M*N*D;
  for(m=0;m<M;m++)for(i=0;i<M;i++) G[m][i]=exp(-dist2(Y[m],Y[i],D)/(2*bet2));

  for(lp=0;lp<nloop;lp++){pren=noise;noise=(pow(2.0*M_PI*sgm2,0.5*D)*M*omg)/(N*(1-omg));
    if(MEM)for(m=0;m<M;m++)for(d=0;d<D;d++) MEM[m+d*M+lp*M*D]=T[m][d];

    /* compute P */
    for(n=0;n<=N;n++) P[M][n]=0;
    for(m=0;m<=M;m++) P[m][N]=0;
    for(m=0;m< M;m++)for(n=0;n<N;n++) P[m][n]=exp(-dist2(X[n],T[m],D)/(2.0*sgm2))+reg;
    for(m=0;m< M;m++)for(n=0;n<N;n++) P[M][n]+=P[m][n];
    for(m=0;m<=M;m++)for(n=0;n<N;n++) P[m][n]/=P[M][n]+noise;
    for(m=0;m< M;m++)for(n=0;n<N;n++) P[m][N]+=P[m][n];
    for(m=0;m< M;m++) P[M][N]+=P[m][N];

    /* compute A and B s.t. AW=B where C=PX */
    for(m=0;m<M;m++)for(d=0;d<D;d++){C[m][d]=0;for(n=0;n<N;n++)C[m][d]+=P[m][n]*X[n][d];}
    for(m=0;m<M;m++)for(i=0;i<M;i++) A[i+m*M]=G[m][i]+(m==i?(lmd*sgm2)/P[m][N]:0);
    for(m=0;m<M;m++)for(d=0;d<D;d++) B[m+d*M]=C[m][d]/P[m][N]-Y[m][d];

    /* solve AW=B and compute T=Y+GW */
    dposv_(&uplo,&M,&D,A,&M,B,&M,&info);
    for(m=0;m<M;m++)for(d=0;d<D;d++) W[m][d]=B[m+d*M];
    for(m=0;m<M;m++)for(d=0;d<D;d++){T[m][d]=Y[m][d];for(i=0;i<M;i++)T[m][d]+=G[m][i]*W[i][d];}

    /* Compute sgm2 */
    pres2=pres1;pres1=sgm2;sgm2=val=0;
    for(n=0;n<N;n++)for(d=0;d<D;d++) val+=SQ(X[n][d])*P[M][n]; sgm2+=val*1;val=0;
    for(m=0;m<M;m++)for(d=0;d<D;d++) val+=C[m][d]*T[m][d];     sgm2-=val*2;val=0;
    for(m=0;m<M;m++)for(d=0;d<D;d++) val+=SQ(T[m][d])*P[m][N]; sgm2+=val*1;val=0;
    sgm2/=P[M][N]*D;

    conv=log(pres2)-log(sgm2 );
    printf("loop=%d\tsgm2=%lf\tnoise=%lf\tconv=%lf\n",lp,sgm2,noise,conv);

    if(fabs(conv)<1e-8)break;
  } NUM=lp==nloop?lp:lp+1;

  return 0;
}

int readPrms(double prms[6],const char *file){
  FILE* fp=fopen(file,"r");
  fscanf(fp,"omg:%lf\n", prms+0);
  fscanf(fp,"bet:%lf\n", prms+1);
  fscanf(fp,"lmd:%lf\n", prms+2);
  fscanf(fp,"nlp:%lf\n", prms+3);
  fscanf(fp,"dz:%lf\n",  prms+4);
  fscanf(fp,"mem:%lf\n", prms+5);
  fclose(fp);
  return 0;
}

int main(int argc, char **argv){
  int M,N,D,nlp,size[3];double dz,prms[6];
  double **W,**T,**G,**P,**C,*A,*B,**X,**Y;

  if(argc==1){printf("./cpd <X> <Y>\n");exit(1);}

  X=readPoints(&N,&D,argv[1]);
  Y=readPoints(&M,&D,argv[2]);
  size[0]=M;size[1]=N;size[2]=D;

  readPrms(prms,"prms.txt"); nlp=prms[3];dz=prms[4];

  rescalePoints(X,N,D,dz);
  rescalePoints(Y,M,D,dz);
  normPoints(X,N,D);
  normPoints(Y,M,D);

  W = calloc2d(M,D); T=calloc2d(M,D);
  G = calloc2d(M,M); C=calloc2d(M,D); P=calloc2d(M+1,N+1);
  A = calloc(M*M,sizeof(double));
  B = calloc(M*D,sizeof(double));
  MEM=prms[5]>0?malloc(3*sizeof(int)+nlp*M*D*sizeof(double)):NULL;

  #define CD (const double **)
  cpd(W,T,G,P,C,A,B,CD X,CD Y,size,prms);

  writePoints("T1.txt",CD T,M,D);
  writePoints("X1.txt",CD X,N,D);
  writePoints("Y1.txt",CD Y,M,D);
  #undef  CD

  if(MEM){FILE *fp=fopen("ontheway.bin","wb");
    fwrite(&M,  sizeof(int),   1,      fp);
    fwrite(&D,  sizeof(int),   1,      fp);
    fwrite(&NUM,sizeof(int),   1,      fp);
    fwrite(MEM, sizeof(double),NUM*M*D,fp);
    fclose(fp);
  }

  return 0;
}

