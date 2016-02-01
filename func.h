int rigid (double       **  W,        /*  D+1 x  D          | Linear map            */
           double       **  T,        /*  M   x  D          | Moved points          */
           double       **  P,        /*  M+1 x  N+1        | Matching probablity   */
           double       *** C,        /*  4 x max(M,N) x D  | Working memory        */
           double       *   S,        /*  nlp x M x D       | Working wemory (3D)   */
           const double **  X,        /*  N   x  D          | Point set 1 (Data)    */
           const double **  Y,        /*  M   x  D          | Point set 2 (Data)    */
           const int        size[3],  /*  M, N, D           | D must be 2 or 3      */
           const double     prms[2],  /*  parameters: nloop, omg                    */
           const int        verb      /*  flag: verbose                             */
);

int affine(double       **  W,        /*  D+1 x  D          | Linear map            */
           double       **  T,        /*  M   x  D          | Moved points          */
           double       **  P,        /*  M+1 x  N+1        | Matching probablity   */
           double       *** C,        /*  4 x max(M,N) x D  | Working memory        */
           double       *   S,        /*  nlp x M x D       | Working wemory (3D)   */
           const double **  X,        /*  N   x  D          | Point set 1 (Data)    */
           const double **  Y,        /*  M   x  D          | Point set 2 (Data)    */
           const int        size[3],  /*  M, N, D           | D must be 2 or 3      */
           const double     prms[2],  /*  parameters: nloop, omg                    */
           const int        verb      /*  flag: verbose                             */
);

int cpd   (double       **  W,        /*  M   x   D         | Displacement matrix   */
           double       **  T,        /*  M   x   D         | Transformed point set */
           double       **  P,        /*  M+1 x   N+1       | Matching probablity   */
           double       **  G,        /*  M   x   M         | Gram matrix of Y      */
           double       *** U,        /*  2 x K+1 x M       | Working wemory (3D)   */
           double       *** V,        /*  2 x M   x D       | Working wemory (3D)   */
           double       *   A,        /*  M   x   M         | Working wemory (1D)   */
           double       *   B,        /*  M   x   D         | Working wemory (1D)   */
           double       *   S,        /*  nlp x M x D       | Working wemory (1D)   */
           const double **  X,        /*  N   x   D         | Point set 1 (Data)    */
           const double **  Y,        /*  N   x   D         | Point set 2 (Data)    */
           const int        size[3],  /*  M,  N,  D         | D must be 2 or 3      */
           const double     prms[5],  /*  parameters: nlp,omg,bet,lmd,rank          */
           const int        verb      /*  flag: verbose                             */
);
