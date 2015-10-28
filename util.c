// Copyright (c) 2014 Osamu Hirose
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
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include"util.h"

#define MAX_LENGTH 256
#define MAX_NWORDS 10
#define INI_SIZE   1024

double dist2(const double x1[3], const double x2[3], const int D){
  assert(D==2||D==3);
  return D==3 ? SQ(x1[0]-x2[0])+SQ(x1[1]-x2[1])+SQ(x1[2]-x2[2])
              : SQ(x1[0]-x2[0])+SQ(x1[1]-x2[1]);
}

double *** calloc3d (const int L, const int M, const int N){
  int l;
  double ***        a    = calloc   (L, sizeof(double**));
  for (l=0;l<L;l++) a[l] = calloc2d (M, N);
  return a;
}

void prepare_sortbox(sortbox *sb, const double * array, const int size){
  int i; for(i=0;i<size;i++){sb[i].val=array[i];sb[i].idx=i;}
  return;
}

int cmp_sortbox(const void *a, const void *b){
  const sortbox *sa = (const sortbox*) a;
  const sortbox *sb = (const sortbox*) b;
  return sa->val>sb->val ? 1 : sa->val<sb->val ? -1 : 0;
}


double ** calloc2d (const int M, const int N){
  int m;
  double **         a    = calloc (M, sizeof(double*));
  for (m=0;m<M;m++) a[m] = calloc (N, sizeof(double ));
  return a;
}

int    ** calloc2i (const int M, const int N){
  int m;
  int **            a    = calloc (M, sizeof(int *));
  for (m=0;m<M;m++) a[m] = calloc (N, sizeof(int  ));
  return a;
}

short ** calloc2s (const int M, const int N){
  int m;
  short **          a    = calloc (M, sizeof(short*));
  for (m=0;m<M;m++) a[m] = calloc (N, sizeof(short ));
  return a;
}

char ** calloc2c (const int uw, const int ulen){
  int w;
  char **            a    = calloc (uw,  sizeof(char*));
  for (w=0;w<uw;w++) a[w] = calloc (ulen,sizeof(char ));
  return a;
}


double ** readPoints(int *N, int *D, const char *file){
  int n; FILE *fp; double **X;
  fp=fopen(file,"rb");
  fread(N,sizeof(int),1,fp);
  fread(D,sizeof(int),1,fp);
  X=calloc2d(*N,*D);
  for(n=0;n<*N;n++)
    fread(X[n],sizeof(double),*D,fp);
  fclose(fp);
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
  int i; if(D==3||fabs(dz-1.0)>1e-8) for(i=0;i<N;i++) X[i][2]*=dz;
  return 0;
}

