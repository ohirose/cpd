#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"version.h"

void printUsage(void){
  printf("\n");
  printf("  USAGE: ./cpd <mode> <X> <Y> (+ options) \n\n");
  printf("  MODE: At least one of characters r, a, and c must be included in <mode>.             \n");
  printf("  Optionally, m and q which specify print options can be attatched.                    \n");
  printf("  r: rigid, a: affine, c: cpd, q: quiet, m: memorize optimization process.           \n\n");
  printf("  X: point set 1, reference points.                                                    \n");
  printf("  Y: point set 2, floating points.                                                   \n\n");
  printf("  OPTIONs: Options must be added AFTER the arguments. If the parameter file is set     \n");
  printf("  as the argument of '-p', other parameters specified by options are ignored.          \n");
  printf("  -n: nloop, -w omega, -l lambda, -b beta, -r rank, -z zscale, -p <parameter file>.  \n\n");
  printf("  EXAMPLE: ./cpd rac X.txt Y.txt -w 0.5 -l 2.0 -b 0.9 -z 3.5 -n 2000 -r 20           \n\n");
}

void printVersion(void){
  printf("\n");
  printf(" cpd version %s (%s).\n",_VERSION_,_DATE_);
  printf(" Copyright (c) Osamu Hirose                                                          \n\n");
  printf(" This software is an implementation of the point-set registration algorithm known as   \n");
  printf(" Coherent Point Drift (CPD) invented by Andriy Myronenko and Xubo Song (2010).         \n");
  printf(" Algorithm details are available in their article \"Point Set Registration: Coherent   \n");
  printf(" Point Drift, IEEE TPAMI, 32(12), 2262--2275, 2010.                                  \n\n");
}

void printInfo(const char **files, const int size[3], const double prms[6], const int flag){
    const char *fX  =files[0]; const int M=size[0];
    const char *fY  =files[1]; const int N=size[1];
    const char *fprm=files[2]; const int D=size[2];
    const char *fout=files[3]; const int K=prms[4];
    printf("\n");
    printf("  Input Data:\n");
    printf("    Point set 1 (reference): [%s]\n", fX);
    printf("    Point set 2 (floating):  [%s]\n", fY);
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

void printOptIndex(const char mode, const int loop, const double Np, const double sigma, const double noise, const double conv){
  char *r="Rigid",*a="Affine",*c="CPD";
  if(loop) printf("\033[F\033[J\033[F\033[J\033[F\033[J\033[F\033[J\033[F\033[J\033[F\033[J");
  printf("  %s: %d loops\n",mode=='r'?r:(mode=='a'?a:c),loop);
  printf("    Np     = %lf\n",  Np);
  printf("    sigma  = %lf\n",  sigma);
  printf("    noise  = %lf\n",  noise);
  printf("    conv   = %lf\n\n",conv);
  return;
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


