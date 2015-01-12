#include <stdio.h>
#include <math.h>
#include "ER.h"

#define pi 3.141592653589793
#define max(a, b) (a >= b? a: b)
#define MIN0(a, b) (a <= b? a: b)
#define COS(a) cos(a)
#define SIN(a) sin(a)


int ZERO(int LX, double *X){
    int I; 
    if (LX <= 0){
        return 1;
    }
    for (I=0; I<LX; I++){ 
        X[I] = 0; 
    }
    return 0; 
}


double DOT(int L, double *X, double *Y){
    /*
    DOT: DOT product
    
    p. 20 
    
    
    P: the dot product 
        
    */
    // check original version 
    double P; 
    int I;
    P = 0.0; 
    if (L <= 0){ 
        return P; 
    }
    for (I=0; I<L; I++){ 
        P += X[I] * Y[I]; 
    }
    return P; 
}


int REMAV(int LY, double *Y){ 
    /*
    REMAV: REMove AVerage
    
    p.22
   */
    double S = 0.;
    int I; 
    for (I=0; I<LY; I++){ 
        S += Y[I]; 
    }
    double AVERAG; 
    AVERAG = S / LY; 
    for (I=0; I<LY; I++){ 
        Y[I] -= AVERAG; 
    }
    return 0; 
}

int CROSS(int LX, double *X, int LY, double *Y, int LG, double *G){
    /*
    CROSS: CROSS correlation
    
    p. 27
    
    */
    int J; 
    ZERO(LG, G); 
    for (J=0; J<LG; J++){
        // print(MIN0((LY, LX - J)))
        // N - J samples
        G[J] = DOT(MIN0(LY, LX - J), &X[J], Y);
//        int i = MIN0(LY, LX - J); 
//        printf("%2i, %3.2f, %3.2f, %3.2f \n", i, X[J], Y[0], G[J]); 
    }
    return 0; 
}

int NORMAG(int LX, double *X){ 
    /*
    NORMAG: NORmalizes an array by the magnitude of the element which is greatest in MAGnitude. 
    
    p.23
    */
    int I; 
    double B; 
    B = 0.0; 
    for (I=0; I<LX; I++){
        // print(B)
        B = max(fabs(X[I]), B); 
    }
    for (I=0; I<LX; I++){ 
        X[I] /= B; 
    }
    return 0; 
}

int MACRO(int N, int LX, double *X, int LY, double *Y, int LG, double *G){ 
    /*
    MACRO: empirical MultichAnnel CROss correlation 

    p. 203
    */
    int I, J, I1, J1, IJ; 
    ZERO(LG * N * N, G); 
    for (I=0; I<N; I++){ 
        I1 = I * LX; 
        for (J=0; J<N; J++){ 
            J1 = J * LY; 
            IJ = LG * N * J + LG * I; 
//            printf("%5i \n", J1); 
            CROSS(LX, &X[I1], LY, &Y[J1], LG, &G[IJ]); 
        }
    }
    return 0; 
}

int COSTAB(int M, double *TABLE){ 
    /*
    full wavelength COSine TABle
    
    p. 60
    */
    int FM, MM, I;  
    FM = M + M - 2; 
    MM = M + M - 1; 
    for (I=0; I<MM; I++){ 
        TABLE[I] = COS(I * 2.0 * pi / FM); 
    }
    return 0; 
}

int SINTAB(int M, double *TABLE){
    /*
    full wavelength SINe TABle
    
    p. 60
    
    */
    int FM, MM, I; 
    FM = M + M - 2; 
    MM = M + M - 1; 
    for (I=0; I<MM; I++){ 
        TABLE[I] = SIN(I * 2.0 * pi / FM);
    }
    return 0; 
}

int COSP(int N, double *DATA, double *TABLE, int M, int K, double *C){
    /*
    COSP: compute kth value of either the COSine transform or sine transform. 
    
    (COSine decomPosition ?)
    
    p.61
    */
    int J, KK, MM, MMM, I; 
    J = 0; 
    C[0] = 0.0; 
    KK = K - 1;
    MM = M + M - 1;
    MMM = MM - 1;
    for (I=0; I<N; I++){
//        printf("%2i, %2i, %2i, %2i\n", I, J, KK, MM); 
        C[0] += DATA[I] * TABLE[J];
        J += KK;
        if ((J + 1 - MM) <= 0){ 
            continue;  
        }
        else {
            J -= MMM; 
        }
    }
    return 0;
}    

int MOVE(int LX, double *X, double *Y){
    /*
    MOVE: MOVE one array from one storage location to another
    
    p. 18
    
    */
    int cond, K, I; 
    // cond = (XLOCF(X) - XLOCF(Y));
    cond = (&X[0] - &Y[0]);
    // print(cond)
    if (cond < 0){
        K = LX - 1; 
        for (I=0; I<LX; I++){
            Y[K] = X[K]; 
            K -= 1; 
        }
    }
    if (cond > 0){ 
        for (I=0; I<LX; I++){ 
            Y[I] = X[I];
        }
    }
    if (cond == 0){
        return 1; 
    }
    return 0; 
}

int QUADCO(int L, int N, double *R, double *S, double *SP){ 
    /*
    QUADCO: QUADrature spectra and COspectra from multichannel autocorrelation function
    
    p. 212
    
    */
    int I, J, K, IJK, IKJ, OneJK; 
    double WEIGHT, ODD, EVEN; 
    for (J=0; J<N; J++){ 
        for (K=J; K<N; K++){ 
            for (I=0; I<L; I++){ 
                WEIGHT = ((double) (L - I)) / (double) L; 
//                 printf("%3.2f\n", WEIGHT); 
                IJK = K * L * N + J * L + I;
                IKJ = J * L * N + K * L + I;
                EVEN = R[IJK] + R[IKJ];
                ODD = R[IJK] - R[IKJ];
//                 printf("%2i, %2i, %2i, %2i, %2i, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f \n", I, J, K, IJK, IKJ, R[IJK], R[IKJ], EVEN, ODD, WEIGHT); 
                R[IKJ] = WEIGHT * ODD; 
                R[IJK] = WEIGHT * EVEN; 
            }
            OneJK = K * L * N + J * L; 
            R[OneJK] /= 2.0;
        }
    }
//    print(R)
    for (J=0; J<N; J++){ 
        for (K=0; K<N; K++){ 
            OneJK = K * L * N + J * L;  
            if (K >= J){ 
                COSTAB(L, SP);
//                printf("%3.2f, %3.2f, %3.2f\n", SP[0], SP[1], SP[2]);
            }
            else { 
                SINTAB(L, SP);
//                printf("%3.2f, %3.2f, %3.2f\n", SP[0], SP[1], SP[2]); 
            }
            for (I=0; I<L; I++){ 
                IJK = K * L * N + J * L + I; 
//                printf("%2i, %2i, %3.2f, %3.2f \n", OneJK, IJK, R[OneJK], S[IJK]); 
                COSP(L, &R[OneJK], SP, L, I + 1, &S[IJK]); 
//                printf("%3.2f\n", S[IJK]); 
            }
            R[OneJK] *= 2.0;
        }
    }
    return 0; 
}


int COHERE(int L, int N, double *S, double *C){
    /** 
    COHERE: COHEREnce
    
    p. 215
    
    **/
    int JP, J, K, I, IJK, IKJ, IJJ, IKK, OneJJ; 
    double num, den, CO, PH;
    for (JP = 1 ; JP < N; JP++){ 
        J = JP - 1; 
        for (K = JP; K < N; K++){ 
            for (I = 0; I < L; I++){
                IJK = L * N * K + L * J + I;
                IKJ = L * N * J + L * K + I;
                IJJ = L * N * J + L * J + I;
                IKK = L * N * K + L * K + I;
                num = pow(S[IJK], 2.) + pow(S[IKJ], 2.);
                den = S[IJJ] * S[IKK];
//                printf("%2i, %2i, %2i \n ", I, J, K);
//                printf("%2i, %2i, %2i, %2i \n", IJK, IKJ, IJJ, IKK); 
//                print(S[IJK], S[IKJ], S[IJJ], S[IKK])
                CO = sqrt(fabs(num / den));
                PH = atan2(S[IKJ], S[IJK]);
//                printf("%3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f\n", S[IKJ], S[IJK], num, den, CO, PH); 
//                PH = ATAN2(S[IJK], S[IKJ])
                C[IJK] = CO;
                C[IKJ] = 180.0 * (PH / pi);
            }
        }
    }
    for (J = 0; J < N; J++){
        OneJJ = L * N * J + L * J; 
//        print(OneJJ)
//        print(C)
//        print(S)
//        print(C)
        MOVE(L, &S[OneJJ], &C[OneJJ]);
//        print(C)
//        print(C[OneJJ:])
        NORMAG(L, &C[OneJJ]);
//        print(C)
    }
    return 0; 
}




