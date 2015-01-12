#include <stdio.h>
#include "ER.h"

#define sunspotX {101, 82, 66, 35, 31, 7, 20, 92, 154, 126, 85, 68, 38, 23, 10, 24, 83, 132, 131, 118, 90, 67, 60, 47, 41, 21, 16, 6, 4, 7, 14, 34, 45, 43, 48, 42, 28, 10, 8, 2, 0, 1, 5, 12, 14, 35, 46, 41, 30, 24, 16, 7, 4, 2, 8, 17, 36, 50, 62, 67, 71, 48, 28, 8, 13, 57, 122, 138, 103, 86, 63, 37, 24, 11, 15, 40, 62, 98, 124, 96, 66, 64, 54, 39, 21, 7, 4, 23, 55, 94, 96, 77, 59, 44, 47, 30, 16, 7, 37, 74, \
155, 113, 3, 10, 0, 0, 12, 86, 102, 20, 98, 116, 87, 131, 168, 173, 238, 146, 0, 0, 0, 0, 12, 0, 37, 14, 11, 28, 19, 30, 11, 26, 0, 29, 47, 36, 35, 17, 0, 3, 6, 18, 15, 0, 3, 9, 64, 126, 38, 33, 71, 24, 20, 22, 13, 35, 84, 119, 86, 71, 115, 91, 43, 67, 60, 49, 100, 150, 178, 187, 76, 75, 100, 68, 93, 20, 51, 72, 118, 146, 101, 61, 87, 53, 69, 46, 47, 35, 74, 104, 97, 106, 113, 103, 68, 67, 82, 89, 102, 110, \
66, 62, 66, 197, 63, 0, 121, 0, 113, 27, 107, 50, 122, 127, 152, 216, 171, 70, 141, 69, 160, 92, 70, 46, 96, 78, 110, 79, 85, 113, 59, 86, 199, 53, 81, 81, 156, 27, 81, 107, 152, 99, 177, 48, 70, 158, 22, 43, 102, 111, 90, 86, 119, 82, 79, 111, 60, 118, 206, 122, 134, 131, 84, 100, 99, 99, 69, 67, 26, 106, 108, 155, 40, 75, 99, 86, 127, 201, 76, 64, 31, 138, 163, 98, 70, 155, 97, 82, 90, 122, 70, 96, 111, 42, 97, 91, 64, 81, 162, 137} 

int printArray(int N, int NMAX, double *A){
    int i; 
    if (N <= NMAX){
        for (i=0; i<N; i++){
            printf("%3.2f ", A[i]);
        }
        printf("\n"); 
    }
    else {
        for (i=0; i<NMAX; i++){
            printf("%3.2f ", A[i]);
        }
        printf("... "); 
        printf("%3.2f \n", A[N-1]);
    }
    return 0; 
}

int printBeg(char *s){
    printf("____________________\n"); 
    printf("%s \n", s);
    return 0; 
}

int printEnd(){
    printf("____________________\n"); 
    return 0; 
}


int testZERO(){
    printBeg("testZERO"); 
    double A[10]; 
    int N = 10; 
    ZERO(N, A); 
    printArray(N, N, A); 
    printEnd(); 
    return 0; 
}

int testDOT(){
    printBeg("testDOT"); 
    double A[3] = {1., 2., 3.}; 
    double B[3] = {2., 3., 4.};
    double P; 
    int N = 3; 
    P = DOT(N, A, B); 
    printf("%3.2f \n", P);
    printEnd(); 
    return 0; 
}

int testCROSS(){
    printBeg("testCROSS"); 
    double A[3] = {1., 2., 3.}; 
    double B[3] = {2., 3., 4.};
    int N = 3; 
    int LG = 3; 
    double G[3];
    CROSS(N, A, N, B, LG, G); 
    printArray(LG, LG, G); 
    printEnd(); 
    return 0; 
}

int testMACRO(){
    printBeg("testMACRO"); 
    double X[300] = sunspotX; 
    int N = 3;
    int LX = 100;
    int LG = 10; 
    double G[N*N*LG];
    printArray(100, 3, &X[0]); 
    printArray(100, 3, &X[100]); 
    printArray(100, 3, &X[200]); 
    REMAV(100, &X[0]);
    REMAV(100, &X[100]);
    REMAV(100, &X[200]);
    printArray(100, 3, &X[0]); 
    printArray(100, 3, &X[100]); 
    printArray(100, 3, &X[200]); 
    MACRO(N, LX, X, LX, X, LG, G); 
    printArray(LG, 3, &G[0 * LG]); 
    printArray(LG, 3, &G[1 * LG]); 
    printArray(LG, 3, &G[2 * LG]); 
    printArray(LG, 3, &G[3 * LG]); 
    printArray(LG, 3, &G[4 * LG]); 
    printArray(LG, 3, &G[5 * LG]); 
    printArray(LG, 3, &G[6 * LG]); 
    printArray(LG, 3, &G[7 * LG]); 
    printArray(LG, 3, &G[8 * LG]); 
    printEnd(); 
    return 0; 
}

int testMOVE(){
    printBeg("testMOVE"); 
    double X[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    printArray(10, 10, X); 
    printf("move 3 from 0 to 1\n"); 
    MOVE(3, &X[0], &X[1]);
    printArray(10, 10, X); 
    printf("move 5 from 2 to 1\n"); 
    MOVE(5, &X[2], &X[1]);
    printArray(10, 10, X); 
    printEnd(); 
    return 0; 
}

int testCOSP(){
    printBeg("testCOSP");
    double TABLE[10]; 
    int M = 4; 
    ZERO(10, TABLE); 
    printf("COSTAB\n"); 
    COSTAB(M, TABLE);
    printArray(10, 10, TABLE);
    printf("SINTAB\n");
    SINTAB(M, TABLE);
    printArray(10, 10, TABLE);
    double DATA[3] = {1, 1, 1}; 
    printf("COSP\n");
    double C[1]; 
    int N = 3; 
    COSTAB(M, TABLE);
    for (int i=0; i<M; i++){
        COSP(N, DATA, TABLE, M, i + 1, C);
        printf("%2i %3.2f \n", i, C[0]);
    }
    printEnd(); 
    return 0; 
}

int testQUADCO(){
    printBeg("testQUADCO");
    int L = 2; 
    int N = 2; 
    double R[8] = {7, -2, -1, 1, -1, -2, 3, 1};
    double S[8]; 
    double SP[2*L-1]; 
    QUADCO(L, N, R, S, SP); 
    printArray(8, 8, S);
    printEnd(); 
    return 0; 
}

int testQUADCO_sunspot(){
    printBeg("testQUADCO_sunspot"); 
    double X[300] = sunspotX; 
    int N = 3;
    int LX = 100;
    int LG = 10; 
    double G[N*N*LG], S[N*N*LG], SP[2*LG-1];
    REMAV(100, &X[0]);
    REMAV(100, &X[100]);
    REMAV(100, &X[200]);
    MACRO(N, LX, X, LX, X, LG, G); 
    QUADCO(LG, N, G, S, SP); 
    for (int i=0; i<N*N; i++){
        printArray(LG, 2, &S[i * LG]);
    }
    printEnd(); 
    return 0; 
}

int testCOHERE(){
    printBeg("testCOHERE");
    int N = 2; 
    int L = 2; 
    double S[8] = {5, 9, 0, -0.0, -1.5, -0.5, 4, 2};
    double C[L*N*N]; 
    COHERE(L, N, S, C); 
    printArray(L*N*N, L*N*N, C); 
    printEnd(); 
    return 0; 
}
int testCOHERE_sunspot(){
    printBeg("testCOHERE_sunspot");
    double X[300] = sunspotX; 
    int N = 3;
    int LX = 100;
    int LG = 10; 
    double G[N*N*LG], S[N*N*LG], SP[2*LG-1];
    REMAV(100, &X[0]);
    REMAV(100, &X[100]);
    REMAV(100, &X[200]);
    MACRO(N, LX, X, LX, X, LG, G); 
    QUADCO(LG, N, G, S, SP);
    COHERE(LG, N, S, S);
    for (int i=0; i<N*N; i++){
        printArray(LG, 2, &S[i * LG]);
    }
    printEnd(); 
    return 0; 
}
int main(){
    testZERO(); 
    testDOT();
    testCROSS(); 
    testMACRO(); 
    testMOVE(); 
    testCOSP(); 
    testQUADCO(); 
    testQUADCO_sunspot();
    testCOHERE(); 
    testCOHERE_sunspot(); 
    return 0; 
}
