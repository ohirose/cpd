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

#include<stdio.h>
#include<assert.h>
#include<math.h>
#include"util.h"
#include"lpk.h"
#define WSIZE 100

double det(const double *A, const int D){
  assert(D==2||D==3);
  return D==2? A[0]*A[3]-A[1]*A[2]:
               A[0+D*0]*(A[1+D*1]*A[2+D*2]-A[1+D*2]*A[2+D*1])
              -A[0+D*1]*(A[1+D*0]*A[2+D*2]-A[1+D*2]*A[2+D*0])
              +A[0+D*2]*(A[1+D*0]*A[2+D*1]-A[1+D*1]*A[2+D*0]);
}

int rot(double       **  W,        /*  D+1 x  D          | Linear map            */
        double       **  T,        /*  M   x  D          | Moved points          */
        double       **  P,        /*  M+1 x  N+1        | Assighment probablity */
        double       *** C,        /*  4 x max(M,N) x D  | Working memory        */
        const double **  X,        /*  N   x  D          | Point set 1 (Data)    */
        const double **  Y,        /*  M   x  D          | Point set 2 (Data)    */
        const int        size[3],  /*  M, N, D           |                       */
        const double     prms[2]   /*  nloop,omg         |                       */
       ){
  int i,j,m,n,d,M,N,D,lp,ws,wi[WSIZE]; int info; char jobz='A';
  int nloop=(int)prms[0]; double omg=prms[1],reg=1e-9;
  double conv,s,noise,sgm2=0,pres1=1e10,pres2=1e20,pren,val,c1,c2;
  double mX[3],mY[3],A[9],B[9],a[3],U[9],S[3],Vt[9],wd[WSIZE];
  double **R,*b,**PXc,**Xc,**Yc,**C1;

  M=size[0];N=size[1];D=size[2]; assert(D<=3);
  R=W;b=W[3];PXc=C[0];Xc=C[1];Yc=C[2];C1=C[3];ws=WSIZE;

  /* initialize */
  for(d=0;d<D;d++)for(i=0;i<D;i++) R[d][i]=d==i?1:0;
  for(m=0;m<M;m++)for(d=0;d<D;d++) T[m][d]=Y[m][d];
  for(m=0;m<M;m++)for(n=0;n<N;n++) sgm2+=dist2(X[n],Y[m],D);sgm2/=M*N*D;

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

    /* centerize X and Y */
    for(d=0;d<D;d++){mX[d]=0;for(n=0;n<N;n++) mX[d]+=X[n][d]*P[M][n];mX[d]/=P[M][N];}
    for(d=0;d<D;d++){mY[d]=0;for(m=0;m<M;m++) mY[d]+=Y[m][d]*P[m][N];mY[d]/=P[M][N];}
    for(n=0;n<N;n++)for(d=0;d<D;d++) Xc[n][d]=X[n][d]-mX[d];
    for(m=0;m<M;m++)for(d=0;d<D;d++) Yc[m][d]=Y[m][d]-mY[d];

    /* compute A and B s.t. AW=B  */
    for(m=0;m<M;m++)for(d=0;d<D;d++){PXc[m][d]=0;for(n=0;n<N;n++) PXc[m][d]+=P  [m][n]*Xc[n][d];}
    for(d=0;d<D;d++)for(i=0;i<D;i++){A[d+i*D ]=0;for(m=0;m<M;m++) A[d+i*D] +=PXc[m][d]*Yc[m][i];}
    for(d=0;d<D;d++)for(i=0;i<D;i++) B[d+i*D]=A[d+i*D];

    dgesdd_(&jobz,&D,&D,B,&D,S,U,&D,Vt,&D,wd,&ws,wi,&info);
    val=det(U,D)*det(Vt,D); if(val<0)for(d=0;d<D;d++)U[d+D*(D-1)]*=-1;
    for(i=0;i<D;i++)for(j=0;j<D;j++){R[i][j]=0;for(d=0;d<D;d++) R[i][j]+=U[i+d*D]*Vt[d+j*D];}

    c1=c2=0;
    for(d=0;d<D;d++)for(i=0;i<D;i++)c1+=A[d+i*D]*R[d][i];
    for(d=0;d<D;d++)for(m=0;m<M;m++)c2+=SQ(Yc[m][d])*P[m][N]; s=c1/c2;

    for(d=0;d<D;d++){val=0;for(i=0;i<D;i++)val+=R[d][i]*mY[i];a[d]=mX[d]-s*val;}
    for(m=0;m<M;m++)for(d=0;d<D;d++){val=0;for(i=0;i<D;i++)val+=R[d][i]*Y[m][i];T[m][d]=s*val+a[d];}

    /* Compute sgm2 */
    pres2=pres1;pres1=sgm2;sgm2=-s*c1;
    for(n=0;n<N;n++)for(d=0;d<D;d++) sgm2+=SQ(X[n][d])*P[M][n]; sgm2/=P[M][N]*D;

    conv=log(pres2)-log(sgm2 );
    printf("loop=%d\tsgm2=%lf\tnoise=%lf\tconv=%lf\n",lp,sgm2,noise,conv);
    if(fabs(conv)<1e-8)break;

  } NUM=lp==nloop?lp:lp+1;

  return 0;
}

