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

int affine(double       **  W,        /*  D+1 x  D          | Linear map            */
           double       **  T,        /*  M   x  D          | Moved points          */
           double       **  P,        /*  M+1 x  N+1        | Assighment probablity */
           double       *** C,        /*  4 x max(M,N) x D  | Working memory        */
           double       *   Q,        /*  nlp x M x D       | Working wemory (3D)   */
           const double **  X,        /*  N   x  D          | Point set 1 (Data)    */
           const double **  Y,        /*  M   x  D          | Point set 2 (Data)    */
           const int        size[3],  /*  M,  N, D          |                       */
           const double     prms[2],  /*  parameters: nloop, omg                    */
           const int        verb      /*  flag: verbose                             */
          ){
  int i,m,n,d,M,N,D,lp,info; char uplo='U';
  int nlp=(int)prms[0]; double omg=prms[1],reg=1e-9;
  double sgm2=0,noise,pres1=1e10,pres2=1e20,conv;
  double mX[3],mY[3],A[9],B[9],a[3];
  double **F,**PXc,**Xc,**Yc,**C1;

  M=size[0];N=size[1];D=size[2]; assert(D<=3);
  F=W;PXc=C[0];Xc=C[1];Yc=C[2];C1=C[3];

  /* initialize */
  for(m=0;m<M;m++)for(d=0;d<D;d++) T[m][d]=Y[m][d];
  for(m=0;m<M;m++)for(n=0;n<N;n++) sgm2+=dist2(X[n],Y[m],D);sgm2/=M*N*D;

  for(lp=0;lp<nlp;lp++){noise=(pow(2.0*M_PI*sgm2,0.5*D)*M*omg)/(N*(1-omg));
    if(Q)for(m=0;m<M;m++)for(d=0;d<D;d++) Q[m+d*M+lp*M*D]=T[m][d];

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
    for(m=0;m<M;m++)for(d=0;d<D;d++){PXc[m][d]=0;for(n=0;n<N;n++) PXc[m][d]+=P[m][n]*Xc[n][d];}
    for(d=0;d<D;d++)for(i=0;i<D;i++) A[d+i*D]=B[d+i*D]=0;
    for(d=0;d<D;d++)for(i=0;i<D;i++)for(m=0;m<M;m++) A[d+i*D]+=Yc [m][d]*Yc[m][i]*P[m][N];
    for(d=0;d<D;d++)for(i=0;i<D;i++)for(m=0;m<M;m++) B[i+d*D]+=PXc[m][d]*Yc[m][i]; //transpose

    /* solve AW=B and compute T=Y+GW */
    dposv_(&uplo,&D,&D,A,&D,B,&D,&info);
    for(d=0;d<D;d++)for(i=0;i<D;i++) F[d][i]=B[i+d*D]; // transpose
    for(d=0;d<D;d++){a[d]=mX[d];for(i=0;i<D;i++) a[d]-=F[d][i]*mY[i];}
    for(m=0;m<M;m++)for(d=0;d<D;d++){T[m][d]=a[d];for(i=0;i<D;i++)T[m][d]+=Y[m][i]*F[d][i];}

    /* Compute sgm2 */
    pres2=pres1;pres1=sgm2;sgm2=0;
    for(m=0;m<M;m++)for(d=0;d<D;d++){C1[m][d]=0;for(i=0;i<D;i++) C1[m][d]+=Yc[m][i]*F[d][i];}
    for(n=0;n<N;n++)for(d=0;d<D;d++) sgm2+=SQ(X[n][d])*P[M][n];
    for(m=0;m<M;m++)for(d=0;d<D;d++) sgm2-=PXc [m][d]*C1[m][d];
    sgm2/=P[M][N]*D;

    conv=log(pres2)-log(sgm2 );
    if(verb) printf("loop=%d\tsgm2=%lf\tnoise=%lf\tconv=%lf\n",lp,sgm2,noise,conv);
    if(fabs(conv)<1e-8)break;
  }

  return lp;
}

