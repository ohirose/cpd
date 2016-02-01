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
#include"lapack.h"
#include"info.h"

int cpd(double       **  W,        /*  M   x   D     | Displacement matrix      */
        double       **  T,        /*  M   x   D     | Transformed point set    */
        double       **  P,        /*  M+1 x   N+1   | Matching probablity      */
        double       **  G,        /*  M   x   M     | Gram matrix of Y         */
        double       *** U,        /*  2 x K+1 x M   | Working wemory (3D)      */
        double       *** V,        /*  2 x M   x D   | Working wemory (3D)      */
        double       *   A,        /*  M   x   M     | Working wemory (1D)      */
        double       *   B,        /*  M   x   D     | Working wemory (1D)      */
        double       *   S,        /*  nlp x M x D   | Working wemory (1D)      */
        const double **  X,        /*  N   x   D     | Point set 1 (Data)       */
        const double **  Y,        /*  N   x   D     | Point set 2 (Data)       */
        const int        size[3],  /*  M,  N,  D     | D must be 2 or 3         */
        const double     prms[5],  /*  parameters: nlp,omg,bet,lmd,rank         */
        const int        verb      /*  flag: verbose                            */
       ){

  int    i,k,m,n,d,M,N,D,K,lp,nlp,info,lwork; char uplo='U',jobz='V';
  double omg,lmd,bet2,reg=1e-9; double **Q,**C0,**C1,**PX,*L,*Lr;
  double sgm2=0,noise,val,pres1=1e10,pres2=1e20,conv;

  nlp=(int)prms[0];omg=prms[1];bet2=SQ(prms[2]);lmd=prms[3];K=prms[4];
  M=size[0];N=size[1];D=size[2];Q=U[0];C0=U[1];PX=V[0];C1=V[1];L=U[0][K];Lr=U[1][K];lwork=M*M;

  /* initialization */
  for(m=0;m<M;m++)for(d=0;d<D;d++) W[m][d]=0;
  for(m=0;m<M;m++)for(d=0;d<D;d++) T[m][d]=Y[m][d];
  for(m=0;m<M;m++)for(n=0;n<N;n++) sgm2+=dist2(X[n],Y[m],D);sgm2/=M*N*D;
  for(m=0;m<M;m++)for(i=0;i<M;i++) G[m][i]=A[i+m*M]=exp(-dist2(Y[m],Y[i],D)/(2*bet2));

  /* additional initialization when low-rank approximation */
  if(K){dsyev_(&jobz,&uplo,&M,A,&M,Lr,B,&lwork,&info);
    for(k=0;k<K;k++) L[k]=Lr[M-k-1];
    for(k=0;k<K;k++)for(m=0;m<M;m++) Q[k][m]=A[m+(M-k-1)*M];
    for(m=0;m<M;m++)for(i=0;i<M;i++){G[m][i]=0;for(k=0;k<K;k++) G[m][i]+=Q[k][m]*Q[k][i]*L[k];}
  }

  /* main computation */
  for(lp=0;lp<nlp;lp++){noise=(pow(2.0*M_PI*sgm2,0.5*D)*M*omg)/(N*(1-omg));
    if(S)for(m=0;m<M;m++)for(d=0;d<D;d++) S[m+d*M+lp*M*D]=T[m][d];

    /* compute matching probability P */
    for(n=0;n<=N;n++) P[M][n]=0;
    for(m=0;m<=M;m++) P[m][N]=0;
    for(m=0;m< M;m++)for(n=0;n<N;n++) P[m][n]=exp(-dist2(X[n],T[m],D)/(2.0*sgm2))+reg;
    for(m=0;m< M;m++)for(n=0;n<N;n++) P[M][n]+=P[m][n];
    for(m=0;m<=M;m++)for(n=0;n<N;n++) P[m][n]/=P[M][n]+noise;
    for(m=0;m< M;m++)for(n=0;n<N;n++) P[m][N]+=P[m][n];
    for(m=0;m< M;m++) P[M][N]+=P[m][N];

    if(K){/* CASE: low-rank, NOTE: C0=d(P1)*Q, C1=id(P1)PX-Y */
      val=lmd*sgm2;
      for(m=0;m<M;m++)for(d=0;d<D;d++){PX[m][d]=0;for(n=0;n<N;n++) PX[m][d]+=P[m][n]*X[n][d];}
      for(m=0;m<M;m++)for(k=0;k<K;k++) C0[k][m]=Q [k][m]*P[m][N];
      for(m=0;m<M;m++)for(d=0;d<D;d++) C1[m][d]=PX[m][d]/P[m][N]-Y[m][d];
      for(k=0;k<K;k++)for(d=0;d<D;d++){B[k+d*K]=0;for(m=0;m<M;m++) B[k+d*K]+=C0[k][m]*C1[m][d];}
      for(k=0;k<K;k++)for(i=0;i<K;i++){A[k+i*K]=0;for(m=0;m<M;m++) A[k+i*K]+=Q [k][m]*C0[i][m];}
      for(k=0;k<K;k++) A[k+k*K]+=val/L[k];
      dposv_(&uplo,&K,&D,A,&K,B,&K,&info); assert(info==0);
      for(m=0;m<M;m++)for(d=0;d<D;d++) W[m][d]=P[m][N]*C1[m][d];
      for(m=0;m<M;m++)for(d=0;d<D;d++)for(k=0;k<K;k++) W[m][d]-=C0[k][m]*B[k+d*K];
      for(m=0;m<M;m++)for(d=0;d<D;d++) W[m][d]/=val;
      for(m=0;m<M;m++)for(d=0;d<D;d++){T[m][d]=Y[m][d];for(i=0;i<M;i++)T[m][d]+=G[m][i]*W[i][d];}
    }
    else{/* CASE: full-rank */
      for(m=0;m<M;m++)for(d=0;d<D;d++){PX[m][d]=0;for(n=0;n<N;n++)PX[m][d]+=P[m][n]*X[n][d];}
      for(m=0;m<M;m++)for(i=0;i<M;i++) A[i+m*M]=G[m][i]+(m==i?(lmd*sgm2)/P[m][N]:0);
      for(m=0;m<M;m++)for(d=0;d<D;d++) B[m+d*M]=PX[m][d]/P[m][N]-Y[m][d];
      dposv_(&uplo,&M,&D,A,&M,B,&M,&info);
      for(m=0;m<M;m++)for(d=0;d<D;d++) W[m][d]=B[m+d*M];
      for(m=0;m<M;m++)for(d=0;d<D;d++){T[m][d]=Y[m][d];for(i=0;i<M;i++)T[m][d]+=G[m][i]*W[i][d];}
    }

    /* compute sgm2 (corresponds to residual) */
    pres2=pres1;pres1=sgm2;sgm2=val=0;
    for(n=0;n<N;n++)for(d=0;d<D;d++) val+=SQ(X[n][d])*P[M][n]; sgm2+=val*1;val=0;
    for(m=0;m<M;m++)for(d=0;d<D;d++) val+=PX[m][d]*T[m][d];    sgm2-=val*2;val=0;
    for(m=0;m<M;m++)for(d=0;d<D;d++) val+=SQ(T[m][d])*P[m][N]; sgm2+=val*1;val=0;
    sgm2/=P[M][N]*D;

    /* check convergence */
    conv=log(pres2)-log(sgm2 );
    if(verb) printOptIndex('c',lp,P[M][N],sqrt(sgm2),noise,conv);
    if(fabs(conv)<1e-8)break;
  }

  return lp;
}
