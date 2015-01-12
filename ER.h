int ZERO(int LX, double *X);
double DOT(int L, double *X, double *Y);
int CROSS(int LX, double *X, int LY, double *Y, int LG, double *G);
int NORMAG(int LX, double *X);
int MACRO(int N, int LX, double *X, int LY, double *Y, int LG, double *G); 
int MOVE(int LX, double *X, double *Y);
int COHERE(int L, int N, double *S, double *C);
int REMAV(int LY, double *Y); 
int COSTAB(int M, double *TABLE);
int SINTAB(int M, double *TABLE);
int COSP(int N, double *DATA, double *TABLE, int M, int K, double *C); 
int QUADCO(int L, int N, double *R, double *S, double *SP); 