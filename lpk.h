
int dsyev_ (char*jobz,char*uplo,int*n,double*a,int*lda,double*w,double*work,int*lwork,int*info);
int dposv_ (char *uplo, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);
int dgesdd_(char *jobz, int *M,int *N,double* A, int* lda, double *S, double *U, int*ldu,
            double* Vt, int *ldvt, double *work, int *lwork, int *iwork, int *info);

