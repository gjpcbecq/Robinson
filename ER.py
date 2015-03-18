"""
Enders A. Robinson
"""
#_______________________________________________________________________________
#_______________________________________________________________________________
import numpy
import math
import matplotlib.pyplot as plt
from temp import *
#_______________________________________________________________________________
#_______________________________________________________________________________
pi = numpy.pi
atan = math.atan
NaN = numpy.nan
sign = numpy.sign
#_______________________________________________________________________________
#_______________________________________________________________________________
PI = pi
XLOCF = id
MIN0 = min
# ATAN2 = math.atan2
COS = math.cos
SIN = math.sin
def MOD(a, b): 
    (q, r) = divmod(a, b)
    return r 
ABS = numpy.abs
CABS = numpy.abs
REAL = numpy.real
AIMAG = numpy.imag
CSQRT = numpy.sqrt
SQRT = numpy.sqrt
ATAN = math.atan
def SQRT(a): 
    return pow(a, 0.5)

#_______________________________________________________________________________
#_______________________________________________________________________________
sunspot_X = numpy.array([
101, 82, 66, 35, 31, 7, 20, 92, 154, 126, 85, 68, 38, 23, 10, 24, 83, 132, 131, 118, 90, 67, 60, 47, 41, 21, 16, 6, 4, 7, 14, 34, 45, 43, 48, 42, 28, 10, 8, 2, 0, 1, 5, 12, 14, 35, 46, 41, 30, 24, 16, 7, 4, 2, 8, 17, 36, 50, 62, 67, 71, 48, 28, 8, 13, 57, 122, 138, 103, 86, 63, 37, 24, 11, 15, 40, 62, 98, 124, 96, 66, 64, 54, 39, 21, 7, 4, 23, 55, 94, 96, 77, 59, 44, 47, 30, 16, 7, 37, 74, 
155, 113, 3, 10, 0, 0, 12, 86, 102, 20, 98, 116, 87, 131, 168, 173, 238, 146, 0, 0, 0, 0, 12, 0, 37, 14, 11, 28, 19, 30, 11, 26, 0, 29, 47, 36, 35, 17, 0, 3, 6, 18, 15, 0, 3, 9, 64, 126, 38, 33, 71, 24, 20, 22, 13, 35, 84, 119, 86, 71, 115, 91, 43, 67, 60, 49, 100, 150, 178, 187, 76, 75, 100, 68, 93, 20, 51, 72, 118, 146, 101, 61, 87, 53, 69, 46, 47, 35, 74, 104, 97, 106, 113, 103, 68, 67, 82, 89, 102, 110, 
66, 62, 66, 197, 63, 0, 121, 0, 113, 27, 107, 50, 122, 127, 152, 216, 171, 70, 141, 69, 160, 92, 70, 46, 96, 78, 110, 79, 85, 113, 59, 86, 199, 53, 81, 81, 156, 27, 81, 107, 152, 99, 177, 48, 70, 158, 22, 43, 102, 111, 90, 86, 119, 82, 79, 111, 60, 118, 206, 122, 134, 131, 84, 100, 99, 99, 69, 67, 26, 106, 108, 155, 40, 75, 99, 86, 127, 201, 76, 64, 31, 138, 163, 98, 70, 155, 97, 82, 90, 122, 70, 96, 111, 42, 97, 91, 64, 81, 162, 137], dtype="float")
#_______________________________________________________________________________
#_______________________________________________________________________________
def AUGURY(N, LX, X, LR, R, SPIKE, FLOOR, LF, F, LY, Y, ERROR):
    """
    AUGURY Multichannel Wiener prediction 
    
    p. 99
    """
    NMAX = 9
    VF = numpy.empty((NMAX * NMAX,))
    VB = numpy.empty((NMAX * NMAX, ))
    R = HEAT(N, 1, LX, X, N, 1, LX, X, LR, R)
    RT = 0.0
    for I in range(N): 
        J = I * N + I
        R[J] *= (1. + SPIKE)
        RT += R[J]
#        print(RT)
    NNLR = N * N * LR
    for L in range(LR): 
        (F, Y, VF, VB) = OMEN(N, L + 1, R, F, Y, VF, VB, Y[NNLR:])
#        print("Y", Y)
        Q = 0.0
        for I in range(N): 
            J = I + I * N
            Q += VF[J]
        ERROR[L] = float(Q) / float(RT)
#        print("ERROR[L], Q, RT : ", ERROR[L], Q, RT)
        LF = L + 1
        if (ERROR[L] <= FLOOR): 
            break
    LY = LX + LF - 1
    Y = BRAINY(N, N, LR, F, N, 1, LX, X, Y)
    return (R, LF, F, LY, Y, ERROR)
#_______________________________________________________________________________
#_______________________________________________________________________________
def OMEN(N, L, R, AF, AB, VF, VB, SP): 
    """
    OMEN 
    
    p. 99
    """
    npr = numpy.round
#    print("AF: ", npr(AF, 3))
#    print("AB: ", npr(AB, 3))
#    print("VF: ", npr(VF, 3))
#    print("VB: ", npr(VB, 3))
#    print("SP: ", npr(SP, 3))
    NMAX = 9
    DF = numpy.empty((NMAX * NMAX, ))
    DB = numpy.empty((NMAX * NMAX, ))
    CF = numpy.empty((NMAX * NMAX, ))
    CB = numpy.empty((NMAX * NMAX, ))
#    print(L)
    if (L == 1):
        VF = MOVE(N * N, R, 0, VF, 0)
        VB = MOVE(N * N, R, 0, VB, 0)
        AF = ZERO(N * N, AF)
        for I in range(N):
            IIZ = I + I * N 
            AF[IIZ] = 1.
#        print(AF)
        AB = MOVE(N * N, AF, 0, AB, 0)
#        print("AF, AB, VF, VB, SP", npr(AF, 3), npr(AB, 3), npr(VF, 3), 
#            npr(VB, 3), npr(SP, 3))
        return (AF, AB, VF, VB)
    ZZO = 1 * N * N
#    print("N, L-1, AB, ZZO, R[ZZO] : ", N, L-1, AB, ZZO, R[ZZO:])
    DB = HEAT(N, N, L - 1, AB, N, N, L - 1, R[ZZO:], 1, DB)
#    print(DB)
    for I in range(N): 
        for J in range(N): 
            IJ = I + J * N
            JI = J + I * N
            DF[IJ] = DB[JI]
    SP = MAINE(N, VB, SP)
#    print(VB, SP)
    CF = BRAINY(N, N, 1, DF, N, N, 1, SP, CF)
#    print("CF : ", CF)
    SP = MAINE(N, VF, SP)
#    print("SP : ", SP)
    CB = BRAINY(N, N, 1, DB, N, N, 1, SP, CB)
    SP = MOVE(N * N * (L - 1), AB, 0, SP, ZZO)
    SP = ZERO(N * N, SP)
    AB = MOVE(N * N * L, SP, 0, AB, 0)
    ZZL = (L - 1) * N * N
    AF[ZZL: ] = ZERO(N * N, AF[ZZL: ])
#    print("AF 1 : ", AF)
#    print("----")
#    print("AF, AB, VF, VB, SP", npr(AF, 3), npr(AB, 3), npr(VF, 3), 
#        npr(VB, 3), npr(SP, 3))
    AB = FORM(N, N, L, AB, CB, AF)
#    print("AB : ", AB)
#    print("AF, CF, SP", AF, CF, SP)
    AF = FORM(N, N, L, AF, CF, SP)
#    print("AF 2 : ", AF)
    VF = FORM(N, N, 1, VF, CF, DB)
#    print("VF : ", VF)
    VB = FORM(N, N, 1, VB, CB, DF)
#    print("VB : ", VB)
#    print("AF, AB, VF, VB, SP", npr(AF, 3), npr(AB, 3), npr(VF, 3), 
#        npr(VB, 3), npr(SP, 3))
#    print("----")
    return (AF, AB, VF, VB)
#_______________________________________________________________________________
#_______________________________________________________________________________
def FORM(M, N, L, A, B, C): 
    """
    FORM form error A(i, j, k) = A(i, j, k) - B(i, ii) * c(ii, j, k)
    
    p. 100
    """
    for I in range(M): 
        for J in range(N): 
            for II in range(N): 
                for K in range(L):
                    IJK = I + J * M + K * M * N
                    III = I + II * M 
                    IIJK = II + J * N + K * N * N
#                    print("A[IJK], B[III], C[IIJK]", A[IJK], B[III], C[IIJK])
                    A[IJK] -= B[III] * C[IIJK]
    return A
#_______________________________________________________________________________
#_______________________________________________________________________________
def ATAN2(Y, X): 
    """    
Description:
ATAN2(Y,X) computes the arctangent of the complex number X + i Y. 
Option:
f95, gnu 
Class:
elemental function 
Syntax:
X = ATAN2(Y,X) 
Arguments:
Y	The type shall be REAL(*). 
X	The type and kind type parameter shall be the same as Y. If Y is zero, then X must be nonzero. 


Return value:
The return value has the same type and kind type parameter as Y. It is the principle value of the complex number X + i Y. If X is nonzero, then it lies in the range -\pi \le \arccos (x) \leq \pi. The sign is positive if Y is positive. If Y is zero, then the return value is zero if X is positive and \pi if X is negative. Finally, if X is zero, then the magnitude of the result is \pi/2. 
Example:
          program test_atan2
            real(4) :: x = 1.e0_4, y = 0.5e0_4
            x = atan2(y,x)
          end program test_atan2
     
    
    ATAN2(0, 0) = not defined
    ATAN2(1, 0) = 0 
    ATAN2(0, 1) = pi / 2 = 1.57
    ATAN2(-1, 0) = pi = 3.14
    ATAN2(0, -1) = -pi / 2 = -1.57
    
    """
    
    if (Y == 0):
        if (X < 0): 
            return -pi
        elif (X == 0):
            return pi / 2
        else : 
            return 0
    else : 
        if (Y > 0): 
            if (X == 0) :
                return pi / 2
            elif (X > 0): 
                return math.atan(Y / X)
            else : 
                return pi + math.atan(Y / X)
        elif (Y < 0): 
            if (X == 0) :
                return - pi / 2
            elif (X > 0) : 
                return math.atan(Y / X)
            else : 
                return - pi + math.atan(Y / X)
    """
    if (y != 0):
        if (x > 0) : 
            z = atan(abs(y / x)) * sign(y)
        elif (x == 0): 
            z = pi / 2. * sign(y)
        else : 
            z = (pi - atan(abs(y / x))) * sign(y)
    else : 
        if (x > 0): 
            z = 0
        elif (x < 0):
            z = pi
        else : # (x == 0) 
            z = NaN
    return z        
    """
#_______________________________________________________________________________
def WIENER_1(N, LX, X, M, LZ, Z, LR, LW, FLOOR, LF, F, E, LY, Y, S): 
    """
    WIENER 
        
    p. 253
    """
    # VERSION 1 OF SUBROUTINE WIENER
    # DIMENSION X(N,LX),Z(M,LX),F(M,N,LR)
    # DIMENSION E(LR),Y(M,LY),S(NN*(5*LR+6)+MN*(LR+2)
    #1+2*M*M
    NN = N * N
    NNLR = NN * LR
    MN = M * R
    # IR = 1
    IR = 0
    # IA = 1 + NNLR
    IA = NNLR
    IB = IA + NNLR
    IAP = IB + NNLR
    IBP = IAP + NNLR
    IVA = IBP + NNLR
    IVB = IVA + NN
    IDA = IVB + NN
    IDB = IDA + NN
    ICA = IDB + NN
    ICB = ICA + NN
    IG = ICB + NN
    ICF = IG + MN * LR
    IGAM = ICF + MN
    IH = IGAM + NN
    IFGT = IH + M * M
    S[IH] = HEAT(M, 1, LZ, Z, M, 1, LZ, Z, 1, S[IH])
    if (LW <= 1): 
        L = LR
    IGZ = IG + MN * LW
    IRZ = IR + NN * LW
    if ((LW >= 1) & (LW < LR)): 
        L = LW
    # if ((LW >= 1) & (LW < LR)) ZERO(MN * (LR - LW), S[IGZ])
        S[IGZ: ] = ZERO(MN * (LR - LW), S[IGZ: ])
    # if ((LW >= 1) & (LW < LR)) ZERO(NN * (LR - LW), S[IRZ])
        S[IRZ: ] = ZERO(NN * (LR - LW), S[IRZ: ])
    # if ((LW >= 1) & (LW >= LR)) L = LR
    if ((LW >= 1) & (LW >= LR)): 
        L = LR
    S[IG: IG + L] = HEAT(M, 1, LZ, Z, N, 1, LX, X, L, S[IG: IG + L])
    S[IR: IR + L] = HEAT(N, 1, LX, X, N, 1, LX, X, L, S[IR: IR + L])
    if ((LW <= 1) | (L <= 1)): 
        # GO TO 2
        pass
    else:
        for K in range(1, L):
            IGK = IG + K * MN
            IRK = IR + K * NN
            # WINDOW = 1.0 - float(K - 1) / float(LW - 1)
            WINDOW = 1.0 - float(K) / float(LW - 1)
            S[IGK: IGK + MN] = SCALE(WINDOW, MN, S[IGK: IGK + MN])
            S[IRK: IRK + NN] = SCALE(WINDOW, NN, S[IRK: IRK + NN])
    # RECUR(N, M, LR, S[IH], S[IG], FLOOR, LF, F, E, S[IA], 
    #1S[IB], S[IAP], S[IBP], S[IVA], S[IVB], S[IDA], S[IDB], S[ICA], S[ICB], 
    #2S[ICF], S[IGAM], S[IFGT])
    RECUR(N, M, LR, S[IH], S[IG], FLOOR, LF, F, E, S[IA], 
        S[IB], S[IAP], S[IBP], S[IVA], S[IVB], S[IDA], S[IDB], S[ICA], S[ICB], 
        S[ICF], S[IGAM], S[IFGT])
    LY = LX + LF - 1
    Y = BRAINY(M, N, LF, F, N, 1, LX, X, Y)
    return (LF, F, E, LY, Y, S) 
#_______________________________________________________________________________
#_______________________________________________________________________________
def RECUR(N, M, LR, R, H, G, FLOOR, LF, F, E, 
    A, B, AP, BP, VA, VB, DA, DB, CA, CB, CF, GAM, FGT): 
    """
    used by WIENER_1
    """
    A = ZERO(N * N * LR, A)
    B = ZERO(N * N * LR, B)
    F = ZERO(N * N * LR, F)
    for I in range(N): 
        for J in range(N): 
            IJ = I + N * J
            IJO = I + N * J 
            VA[IJ] = R[IJO]
        IIO = I + I * N
        A[IIO] = 1.
        B[IIO] = 1.
    G = SIMEQ1(M, N, F, R, G)
    LF = 1
    FGT = HEAT(M, N, 1, F, M, N, 1, G, 1, FGT)
    E[0] = 1.0 - SPUR(M, FGT) / SPUR(M, H)
    if (E[0] <= FLOOR):
        return 
    if (LR == 1): 
        return 
    for L in range(1, LR): 
        DA = ZERO(N * N, DA)
        OOL = L * M * N 
        GAM = MOVE(M * N, G, OOL, GAM, 0)
        for I in range(N): 
            for LI in range(L): 
                # LD = L - LI + 1
                LD = L - LI + 1
                for K in range(N):
                    for J in range(N): 
                        IJ = I + N * J
                        IKLI = I + N * K + N * N * LI
                        DA[IJ] -= A[IKLI] * R[KJLD]
                    for J in range(M): 
                        JI = J + M * I
                        JKLI = J + M * K + M * N * LI
                        KILD = K + M * I + M * N * LD
                        GAM[JI] -= F[JKLI] * R[KILD]
            # 4 end for LI
            for J in range(N):
                IJ = I + N * J
                JI = J + N * I
                DB[JI] = DA[IJ]
        # 5 end for I 
#_______________________________________________________________________________
#_______________________________________________________________________________
def ZERO(LX, X):
    if (LX <= 0): 
        return
    for I in range(LX): 
        X[I] = 0
    return X
#_______________________________________________________________________________
#_______________________________________________________________________________
def MOVE(LX, X, IX, Y, IY):
    """
    MOVE: MOVE one array from one storage location to another
    
    p. 18
    
    """
    # print("LOC X IX, LOC Y IY", XLOCF(X), XLOCF(Y))
    cond = (XLOCF(X) - XLOCF(Y))
    # print(cond)
    if (cond < 0):
        K = LX - 1
        for I in range(0, LX): 
            Y[IY + K] = X[IX + K].copy()
            K -= 1
    elif (cond > 0): 
        for I in range(0, LX): 
            Y[IY + I] = X[IX + I].copy()
    else : # (cond == 0)
        cond = (IX - IY)
        if (cond < 0): 
            K = LX - 1
            for I in range(0, LX): 
                Y[IY + K] = X[IX + K].copy()
                K -= 1
        elif (cond > 0):
            for I in range(0, LX): 
                Y[IY + I] = X[IX + I].copy()
        else : 
            pass
    return Y
#_______________________________________________________________________________
#_______________________________________________________________________________
def DOT(L, X, Y):
    """
    DOT: DOT product
    
    p. 20 
    
    
    P: the dot product 
        
    """
    P = 0.0
    if (L <= 0) : 
        return P 
    for I in range(L): 
        P += X[I] * Y[I]
    return P 
#_______________________________________________________________________________
#_______________________________________________________________________________
def DOTR(L, X, Y):
    """
    DOTR: DOT product reverse of two vectors
    
    p. 20 
    
    P: the dot product 
    """
    P = 0.0
    if (L <= 0) : 
        return P 
    for I in range(L): 
        J = L - I - 1
        P += X[I] * Y[J]
    return P 
#_______________________________________________________________________________
#_______________________________________________________________________________
def REVERS(LX, X): 
    """
    REVERS reverse the order of the element
    
    p.22
    """
    L = LX / 2
    for I in range(L):
        J = LX - I - 1
        TEMP = X[I]
        X[I] = X[J]
        X[J] = TEMP
    return X
#_______________________________________________________________________________
#_______________________________________________________________________________
def CROSS(LX, X, LY, Y, LG, G):
    """
    CROSS: CROSS correlation
    
    p. 27
    
    """
    #G = numpy.zeros((LG, ))
    for J in range(LG):
        # print(MIN0((LY, LX - J)))
        # N - J samples
        # MIN0((LY, LX - (J + 1) + 1))
        G[J] = DOT(MIN0((LY, LX - J)), X[J:], Y)
    return G
#_______________________________________________________________________________
#_______________________________________________________________________________
def REMAV(LY, Y): 
    """
    REMAV: REMove AVerage
    
    p.22
    """
    S = 0.
    for I in range(LY): 
        S += Y[I]
    AVERAG = float(S) / float(LY)
    for I in range(LY): 
        Y[I] -= AVERAG
    return (Y, AVERAG)
#_______________________________________________________________________________
#_______________________________________________________________________________
def NORMAG(LX, X): 
    """
    NORMAG: NORmalizes an array by the magnitude of the element which is greates in MAGnitude. 
    
    p.23
    """
    B = 0.0
    for I in range(LX):
        # print(B)
        B = max(abs(X[I]), B)
    for I in range(LX): 
        X[I] /= B
    return X
#_______________________________________________________________________________
#_______________________________________________________________________________
def MACRO(N, LX, X, LY, Y, LG, G): 
    """
    MACRO: empirical MultichAnnel CROss correlation 

    p. 203
    """
    G = numpy.zeros((LG * N * N, ))
    for I in range(N): 
        I1 = I * LX
        for J in range(N): 
            J1 = J * LY
            IJ = LG * I + LG * N * J 
            G[IJ:IJ+LG] = CROSS(LX, X[I1:], LY, Y[J1:], LG, G[IJ:IJ+LG])
    return G
#_______________________________________________________________________________
#_______________________________________________________________________________
def COSTAB(M): 
    """
    full wavelength COSine TABle
    
    p. 60
    """
    FM = M + M - 2
    MM = M + M - 1
    TABLE = numpy.empty((MM, ))
    for I in range(MM): 
        TABLE[I] = COS(float(I) * 2 * pi / FM)
    return TABLE
#_______________________________________________________________________________
#_______________________________________________________________________________
def SINTAB(M): 
    """
    full wavelength SINe TABle
    
    p. 60
    
    """
    FM = M + M - 2
    MM = M + M - 1
    TABLE = numpy.empty((MM, ))
    for I in range(MM): 
        TABLE[I] = SIN(float(I) * 2 * pi / FM)
    return TABLE
#_______________________________________________________________________________
#_______________________________________________________________________________
def POLYDV(N, DVS, M, DVD, L, Q): 
    """
    POLYDV divide one polynomial by an other, that is, to deconvolve one signal by another. 
    
    p. 31
    """
    Q = numpy.zeros((L, ))
    Q = MOVE(MIN0(M, L), DVD, 0, Q, 0)
    for I in range(L):
        Q[I] /= DVS[0]
        if (I == L - 1): 
            return 
        K = I
        ISUB = MIN0(N - 1, L - I)
        for J in range(ISUB): 
            K += 1
            Q[K] -= Q[I] * DVS[J + 1]
    return Q
#_______________________________________________________________________________
#_______________________________________________________________________________
def PSQRT(N, C, M, A):
    """
    PSQRT finds the first m + 1 coefficients of the square-root power series
    
    $$ (c_0 + c_1 \, z + \cdots + c_n \, z^n) ^{0.5} = a_0 + a_1 \, z +
     + a_2 \, z^2 + \cdots $$
     
    p. 32
    """
    A[0] = SQRT(C[0])
    TA = 2 * A[0]
    A[1] = C[1] / TA
    A[2] = (C[2] - A[1] * A[1]) / TA
    for I in range(3, M): 
        if (I > N): 
            PA = 0.
        else : 
            PA = C[I]
        PS = 0.
        IH = I / 2
        for J in range(1, IH): 
            K = I - J
            PS += A[J] * A[K]
        PA -= 2. * PS
        if (2 * IH != I): 
            PA -= A[IH] * A[IH]
        A[I] = PA / TA
    return A
#_______________________________________________________________________________
#_______________________________________________________________________________
def POLYEV(LA, A, Z, AZ): 
    """
    POLYEV computes the value of the polynomial at the complex point z
    
    p. 33
    """
    AZ = complex(0.0, 0.0)
    for I in range(LA): 
        J = LA - I - 1
        AZ = Z * AZ + A[J]
    return AZ
#_______________________________________________________________________________
#_______________________________________________________________________________
def CAST(LW, W, LT, TCOS, TSIN, AMP, PHZ): 
    """
    CAST Cosine And Sine Transform
    
    p. 67
    """
    A = numpy.empty((LT, ))
    DANG = pi / float(LT - 1)
    ANG = 0.0
    for I in range(LT): 
        X = complex(COS(ANG), SIN(ANG))
        A = POLYEV(LW, W, X, A)
        TCOS[I] = REAL(A)
        TSIN[I] = AIMAG(A)
        ANG += DANG
    (AMP, PHZ) = POLAR(LT, TCOS, TSIN, AMP, PHZ)
    PHZ = DRUM(LT, PHZ)
    return (TCOS, TSIN, AMP, PHZ)
#_______________________________________________________________________________
#_______________________________________________________________________________
def POLAR(L, RE, XIM, AMP, PHZ):
    """
    POLAR computes polar coordinates
    
    p. 66
    """
    for I in range(L): 
        AMP[I] = SQRT(RE[I] ** 2 + XIM[I] ** 2)
        if (XIM[I] < 0): 
            if (RE[I] < 0): 
                PHZ[I] = ATAN(XIM[I] / RE[I]) - PI
            elif (RE[I] == 0): 
                PHZ[I] = - PI / 2.0
            else: 
                PHZ[I] = ATAN(XIM[I]  / RE[I])
        elif (XIM[I] == 0): 
            if (RE[I] < 0): 
                PHZ[I] = - PI
            elif (RE[I] == 0): 
                PHZ[I] = 0.0
            else: 
                PHZ[I] = ATAN(XIM[I]  / RE[I])
        else: 
            if (RE[I] < 0): 
                PHZ[I] = ATAN(XIM[I]  / RE[I]) + PI
            elif (RE[I] == 0): 
                PHZ[I] = PI / 2.0
            else: 
                PHZ[I] = ATAN(XIM[I]  / RE[I])
    return (AMP, PHZ)
#_______________________________________________________________________________
#_______________________________________________________________________________
def DRUM(LPHZ, PHZ): 
    """
    DRUM makes a phase curve continuous
    
    p.65
    """
    PJ = 0
    for I in range(1, LPHZ): 
        test1 = ABS(PHZ[I] + PJ - PHZ[I - 1]) - PI
        if (test1 > 0):
            test2 = (PHZ[I] + PJ - PHZ[I - 1])
            if (test2 < 0): 
                PJ += 2 * PI
                PHZ[I] += PJ
            elif (test2 == 0): 
                PHZ[I] += PJ
            else: 
                PJ -= 2 * PI
                PHZ[I] += PJ
        else: 
            PHZ[I] += PJ
    return PHZ
#_______________________________________________________________________________
#_______________________________________________________________________________
def COSP(N, DATA, TABLE, M, K): 
    """
    COSP: compute kth value of either the COSine transform or sine transform. 
    
    (COSine decomPosition ?)
    
    p.61
    """
    J = 0
    C = 0.0
    KK = K - 1
    MM = M + M - 1
    MMM = MM - 1
    for I in range(N):
#        print((I, J, KK, MM))
        C += DATA[I] * TABLE[J]
        J += KK
        if ((J + 1 - MM) <= 0): 
            continue 
        else : 
            J -= MMM
    return C
#_______________________________________________________________________________
#_______________________________________________________________________________
def COSTR(LR, R, W, S): 
    """
    COSTR one value of cosine transform of two-sided function
    
    p. 90 
    """
    COSNW = 1.
    SINNW = 0.
    COSW = COS(W)
    SINW = SIN(W)
    S = R[0]
    for I in range(1, LR): 
        T = COSW * COSNW - SINW * SINNW
        COSNW = T
        S += 2 * R[I] * COSNW
    return S
#_______________________________________________________________________________
#_______________________________________________________________________________
def EUREKA(LR, R, G, F, A): 
    """
    EUREKA solve single-channel normal equation by recursion
    
    p. 44
    """
    V = R[0]
    D = R[1]
    A[0] = 1.
    F[0] = G[0] / V
    Q = F[0] * R[1]
    if (LR == 1): 
        return (F, A)
    for L in range(1, LR): 
#        print("---------------")
#        print(("L", L))
        A[L] = - D / V
#        if ((L + 1) != 2):
        if (L != 1):
#            L1 = ((L + 1) - 2) / 2
            L1 = (L - 1) / 2
            L2 = L1 + 1
            if (L2 >= 2):
                for J in range(1, L2): 
                    HOLD = A[J]
#                    K = (L + 1) - (J + 1) + 1 - 1
                    K = L - J
                    A[J] += A[L] * A[K]
                    A[K] += A[L] * HOLD
#                    print("<< J, L, K", J, L, K) 
#            if (2 * L1 != (L + 1) - 2): 
            if (2 * L1 != L - 1): 
#                A[(L2 + 1) - 1] += A[L] * A[(L2 + 1) - 1]
                A[L2] += A[L] * A[L2]
#                print("L2, L, L2",L2, L, L2)
        V += A[L] * D
        F[L] = (G[L] - Q) / V
#         L3 = (L + 1) - 1
        L3 = L 
        for J in range(L3): 
#             K = (L + 1) - (J + 1) + 1 - 1
            K = L - J
            F[J] += F[L] * A[K]
#            print(("J, L, K", J, L, K))
#        print((LR, L, numpy.round(F, 4)))
        if (L == LR - 1): 
            return (F, A)
        D = 0.0
        Q = 0.0
        for I in range(L + 1):
#            K = (L + 1) - (I + 1) + 2 - 1
            K = L - I + 1
#            print(("I, K", I, K))
            D += A[I] * R[K]
            Q += F[I] * R[K] 
    return (F, A)
#_______________________________________________________________________________
#_______________________________________________________________________________
def IMPULS(LD, D, K): 
    """
    IMPULS store an impulse function in an array
    
    p. 26
    """
    for I in range(LD): 
        D[I] = 0.0
    D[K] = 1.0
    return D
#_______________________________________________________________________________
#_______________________________________________________________________________
def INVTOP(LR, R, Q, SPACE): 
    """
    INVTOP inverse of a Toplitz matrix
    
    p. 46
    """
    for K in range(LR): 
#        print("K: {0}".format(K))
        SPACE = IMPULS(LR, SPACE, K)
        J = LR * K
        (Q[J: J + LR], SPACE[LR: 2 * LR]) = EUREKA(LR, R, SPACE, 
            Q[J: J + LR], SPACE[LR: 2 * LR])
    return (Q, SPACE)
#_______________________________________________________________________________
#_______________________________________________________________________________
def SHAPE(LB, B, LD, D, LA, A, LC, C, ASE, SPACE):
    """
    SHAPE computes the waveshaping filter
    
    p. 75
    """
    SPACE = CROSS(LB, B, LB, B, LA, SPACE)
    SPACE[LA: ] = CROSS(LD, D, LB, B, LA, SPACE[LA:])
    (A, SPACE[2 * LA: ]) = EUREKA(LA, SPACE, SPACE[LA: 2 * LA], A, 
        SPACE[2 * LA: ])
    DD = DOT(LD, D, D)
    AG = DOT(LA, A, SPACE[LA: ])
    ASE = (DD - AG) / DD
    (LC, C) = FOLD(LA, A, LB, B, LC, C)
    return (A, LC, C, ASE, SPACE)
#_______________________________________________________________________________
#_______________________________________________________________________________
def SHAPER(LB, B, LD, D, LA, A, LC, C, INDEX, ERRORS, S):
    """
    SHAPER computes waveshaping filter 
    
    p. 83
    """
    print("LB", LB, "B", B, "LD", LD, "D", D, "LA", LA, "A", A, 
        "LC", LC, "C", C, "INDEX", INDEX, "ERRORS", ERRORS, "S", S)
    LC = LA + LB - 1
    LCD = LC + LD - 1
    DD = DOT(LD, D, D)
    S = CROSS(LB, B, LB, B, LA, S)
    for I in range(LCD): 
        C = ZERO(LCD, C)
#        LDI + 1 = LD - (I + 1) + 1
        LDI = LD - I - 1
#        if (I + 1 <= LD):
        if (I <= LD - 1):
            C = MOVE(I, D, LDI, C, 0)
#        ILD + 1 = I + 1 - LD + 1
        ILD = I - LD + 1
#        if (I + 1 >= LD): 
        if (I >= LD - 1): 
            C = MOVE(LD, D, 0, C, ILD)
        S[LA: ] = CROSS(LC, C, LB, B, LA, S[LA: ])
#        if (I + 1 >= 2): 
        if (I >= 1): 
            A = SIDE(S[LA], LA, A, S[2 * LA: ], S)
        else: 
            (A, S[2 * LA: ]) = EUREKA(LA, S, S[LA:], A, S[2 * LA: ])
        AG = DOT(LA, A, S[LA: ])
        (LC, C[:LC]) = FOLD(LA, A, LB, B, LC, C)
        ERRORS[I] = (DD - AG) / DD
    (EMIN, INDEX) = MINSN(LCD, ERRORS, 0, INDEX)
    C = ZERO(LCD, C)
#    LDIND + 1 = LD - (INDEX + 1) + 1
    LDIND = LD - INDEX - 1
#    if (INDEX + 1 <= LD): 
    if (INDEX <= LD - 1): 
        C = MOVE(INDEX, D, LDIND, C, 0)
#    INDLD + 1 = (INDEX + 1) - LD + 1
    INDLD = INDEX - LD + 1
#    if (INDEX + 1 >= LD): 
    if (INDEX >= LD - 1): 
        C = MOVE(LD, D, 0, C, INDLD)
    (A, LC, C, EMIN, S) = SHAPE(LB, B, LC, C, LA, A, LC, C, EMIN, S)
    return (A, LC, C, INDEX, ERRORS, S)
#_______________________________________________________________________________
#_______________________________________________________________________________
def SPIKE(LB, B, LA, A, INDEX, ASE, S): 
    """
    SPIKE computes the spiking filter for the optimum spike position
    
    p. 79
    """
    LD = LA + LB - 1
    print(LD)
    for I in range(LD): 
        S[LD: ] = IMPULS(LD, S[LD: ], I)
        (A, LD, S[LD: 2 * LD], S[I], S[2 * LD: ]) = SHAPE(LB, 
            B, LD, S[LD: ], LA, A, LD, S[LD: ], S[I], S[2 * LD: ])
#        print(I, S)
#        print(A)
    (ASE, INDEX) = MINSN(LD, S, ASE, INDEX)
#    print(ASE, INDEX)
    S[LD: ] = IMPULS(LD, S[LD: ], INDEX)
#    print(S)
    (A, LD, S[LD: 2 * LD], ASE, S[2 * LD: ]) = SHAPE(LB, 
        B, LD, S[LD: ], LA, A, LD, S[LD: ], ASE, S[2 * LD: ])
    return (A, INDEX, ASE, S)
#_______________________________________________________________________________
#_______________________________________________________________________________
def SIDE(H, LF, F, A, R): 
    """
    SIDE Simpson sideways recursion 
    
    p. 80
    """
    V = R[0]
    S = 0.0
    T = 0.0
#    print("A", A, "R", R)
    if (LF != 1): 
        for I in range(1, LF): 
#            J + 1 = LF + 2 - (I + 1) 
            J = LF - I
            S += F[I - 1] * R[I]
            T += A[J] * R[I]
            V += A[I] * R[I]
    FLF = F[LF - 1]
    W = (H - S + FLF * T) / V
    if (LF != 1): 
        for I in range(1, LF):
#            J + 1 = LF - (I + 1) + 2 
            J = LF - I
#            print("I", I, "J", J, "F", F, "W", W, "A", A, "FLF", FLF)
            F[J] = F[J - 1] + W * A[J] - FLF * A[I]
    F[0]= W 
    return F
#_______________________________________________________________________________
#_______________________________________________________________________________
def SPIKER(LB, B, LA, A, LC, C, INDEX, ERRORS, SPACE): 
    """
    SPIKER spiking filter with Simpson's sideways recursion
    
    p. 82
    """
#    print(LB, B, LA, A, LC, C, INDEX, ERRORS, SPACE)
    LC = LA + LB - 1
    SPACE = CROSS(LB, B, LB, B, LA, SPACE)
    for I in range(LC): 
        C = IMPULS(LC, C, I)
        SPACE[LA: ] = CROSS(LC, C, LB, B, LA, SPACE[LA: ])
#        if ((I + 1) >= 2)
        if (I >= 1): 
            A = SIDE(SPACE[LA], LA, A, SPACE[2 * LA: ], SPACE)
        else: 
            (A, SPACE[2 * LA: ]) = EUREKA(LA, SPACE, SPACE[LA: ], A, 
                SPACE[2 * LA: ])
        Q = DOT(LA, A, SPACE[LA:])
        (LC, C) = FOLD(LA, A, LB, B, LC, C)
        ERRORS[I] = 1. - Q
    (EMIN, INDEX) = MINSN(LC, ERRORS, 0, 0)
    C = IMPULS(LC, C, INDEX)
    (A, LC, C, EMIN, SPACE) = SHAPE(LB, B, LC, C, LA, A, LC, C, EMIN, SPACE)
    return (A, LC, C, INDEX, ERRORS)
#_______________________________________________________________________________
#_______________________________________________________________________________
def SPUR(N, A):
    """
    SPUR trac or spur of the diagonal elements of a matrix n x n
    
    p. 49
    """
    SP = 0.0
    for I in range(N): 
        SP += A[I, I]
    return SP
#_______________________________________________________________________________
#_______________________________________________________________________________
def MATRAN(M, N, MATRIX): 
    """
    Matrix transposition
    
    p. 52
    """
    K = M * N - 1
    """
    something about the right bits
    this seems to enable to control if the element has already been changed so that it avoid to repeat the permutation again. 
    I put another matrix filled with 1 to check even if I multiply the space by 2. 
    
    for I in range(1, K): 
        MATRIX[I] = (MATRIX[I] / 2) * 2
    for L in range(2, K): 
        if (MATRIX(L) != (MATRIX[I] / 2) * 2): 
            return 
    """
    CTRL = numpy.zeros(M * N)
    CTRL[0] = 1
    CTRL[-1] = 1
    for L in range(1, K):
        if CTRL[L] == 1: 
            continue
        KEEP = MATRIX[L]
        IJ = L
        cond = True
        while cond: 
#            JLESS1 = IJ // M
#           J is JLESS1 in python or C
            J = IJ // M
            I = IJ - J * M
#            J = JLESS1
            JI = J + I * N
#            print("L, IJ, I, J, JI", L, IJ, I, J, JI)
            KATCH = MATRIX[JI]
            MATRIX[JI] = KEEP
#            print("KEEP, KATCH", KEEP, KATCH)
            KEEP = KATCH
#            print("MATRIX", MATRIX)
            CTRL[IJ] = 1
            IJ = JI
            cond = (IJ != L)
#            print(cond)
    return MATRIX 
#_______________________________________________________________________________
#_______________________________________________________________________________
def QUADCO(L, N, R): 
    """
    QUADCO: QUADrature spectra and COspectra from multichannel autocorrelation function
    
    p. 212
    
    """
    S = numpy.empty((L * N * N, ))
    SP = numpy.empty((2 * L - 1, ))
    for J in range(N): 
        for K in range(J, N): 
            for I in range(L): 
                WEIGHT = float(L - I) / float(L)
#                 print(WEIGHT)
                IJK = L * N * K + L * J + I
                IKJ = L * N * J + L * K + I
                EVEN = R[IJK] + R[IKJ]
                ODD = R[IJK] - R[IKJ]
#                 print((I, J, K, IJK, IKJ, R[IJK], R[IKJ], EVEN, ODD, WEIGHT))
                R[IKJ] = WEIGHT * ODD
                R[IJK] = WEIGHT * EVEN
            OneJK = L * N * K + L * J
            R[OneJK] /= 2.0
    # print(R)
    SPCOS = COSTAB(L)
    SPSIN = SINTAB(L)
    for J in range(N): 
        for K in range(N): 
            OneJK = L * N * K + L * J 
            if (K >= J) : 
                SP = SPCOS
#                 print(SP)
            else : 
                SP = SPSIN
#                 print(SP)
            for I in range(L): 
                IJK = L * N * K + L * J + I
#                 print((OneJK, IJK))
                S[IJK] = COSP(L, R[OneJK:], SP, L, I + 1)
            R[OneJK] *= 2.0
    return S
#_______________________________________________________________________________
#_______________________________________________________________________________
def HEAT(NRX, NCX, LX, X, NRY, NCY, LY, Y, LG, G): 
    """
    HEAT multichannel autocorrelation or cross correlation
    
    G = HEAT(NRX, NCX, LX, X, NRY, NCY, LY, Y, LG, G)
    
    p.207
    """
#    ZERO(NRX * NRY * LG, G)
#    print(NRX, NRY, LG)
    G = numpy.zeros((NRX * NRY * LG, ))
    MIN = MIN0(LG, LX)
    for M in range(NRX): 
        for N in range(NRY):
            for L in range(NCX): 
                for J in range(MIN): 
                    LDOT = MIN0(LY, LX - J)
                    for I in range(LDOT): 
#                        print("M, N, L, J, I : ", M, N, L, J, I)
                        K = I + J 
                        MNJ = M + N * NRX + J * NRX * NRY
                        MLK = M + L * NRX + K * NRX * NCX
                        NLI = N + L * NRY + I * NRY * NCX
                        G[MNJ] += X[MLK] * Y[NLI]
    return G
#_______________________________________________________________________________
#_______________________________________________________________________________
def HEAT_CPLX(NRX, NCX, LX, X, NRY, NCY, LY, Y, LG, G): 
    """
    HEAT multichannel autocorrelation or cross correlation
    
    p.207
    """
#    ZERO(NRX * NRY * LG, G)
    X = X.astype("complex")
    Y = Y.astype("complex")
    G = numpy.zeros((NRX * NRY * LG, ), "complex")
    MIN = MIN0(LG, LX)
    for M in range(NRX): 
        for N in range(NRY):
            for L in range(NCX): 
                for J in range(MIN): 
                    LDOT = MIN0(LY, LX - J)
                    for I in range(LDOT): 
                        K = I + J
                        MNJ = M + N * NRX + J * NRX * NRY
                        MLK = M + L * NRX + K * NRX * NCX
                        NLI = N + L * NRY + I * NRY * NCX
                        G[MNJ] += X[MLK] * Y[NLI].conjugate()
    return G
#_______________________________________________________________________________
#_______________________________________________________________________________
def POLRT(XCOF, COF, M, ROOTR, ROOTI, IER): 

    """
    POLRT compute the real and complex roots of a polynomial
    
    p. 35
    
    """
    ROOTR = numpy.zeros((M, ))
    ROOTI = numpy.zeros((M, ))
    XO = 0
    YO = 0
    IN = 0
    TEMP = 0.0
    IFIT = 0
    N = M 
    IER = 0
    if (XCOF[N] == 0):
        # set error code to 4
        IER = 4
        return (COF, ROOTR, ROOTI, IER)
    if (N <= 0): 
        # set error code to 1
        IER = 1
        return (COF, ROOTR, ROOTI, IER)
    if (N > 36): 
        # set error code to 2
        IER = 2
        return (COF, ROOTR, ROOTI, IER)
    NX = N
    NXX = N + 1
    N2 = 0
    KJ1 = N + 1
    ICT = 0
    SUMSQ = 0.
    ALPHA = 0.
    for L in range(KJ1): 
        # (MT + 1) = KJ1 - (L + 1) + 1
        MT = KJ1 - L - 1
        COF[MT] = XCOF[L]
    # print("XCOF, COF", XCOF, COF)
    (XO, YO, IN) = POLRT_set_initial_values(XO, YO, IN)
    (XO, YO, X, Y, ICT, IN) = POLRT_increment_initial_values(XO, YO, ICT, IN)
    while (N > 0): 
        # print((N, X, Y))
        # evaluate polynomial and derivatives 
        (flag, X, Y, SUMSQ, N, DX, DY) = POLRT_eval_pol_and_deriv(COF, X, Y, SUMSQ, N)
        # print(flag)
        if (flag == 0): 
            X = 0.0
            NX -= 1
            NXX -= 1
            Y = 0.0
            SUMSQ = 0.0
            ALPHA = X
            N -= 1
            (COF, ROOTI, ROOTR, N2) = POLRT_compute_root(COF, ALPHA, N, SUMSQ, 
                ROOTI, ROOTR, N2, X, Y)
            (XO, YO, IN) = POLRT_set_initial_values(XO, YO, IN)
            (XO, YO, X, Y, ICT, IN) = POLRT_increment_initial_values(XO, YO, ICT, IN)
        else: 
            if (SUMSQ == 0):
                if (IFIT == 0): 
                    (XO, YO, X, Y, ICT, IN) = POLRT_increment_initial_values(XO, YO, ICT, IN)
                else : 
                    X = XPR
                    Y = YPR
                    IFIT = 0
                    (X, Y, N, SUMSQ, ALPHA) = POLRT_check_root(X, Y, N)
                    (COF, ROOTI, ROOTR, N2) = POLRT_compute_root(COF, ALPHA, N, 
                        SUMSQ, ROOTI, ROOTR, N2, X, Y)
                    (XO, YO, IN) = POLRT_set_initial_values(XO, YO, IN)
                    (XO, YO, X, Y, ICT, IN) = POLRT_increment_initial_values(XO, YO, ICT, IN)
            else: 
                if (ABS(DY) + ABS(DX) >= 1e-5):
                    # step iteration counter
                    ICT += 1
                    # print(ABS(DY) + ABS(DX))
                    if (ICT < 500): 
                        # print(SUMSQ)
                        # print(ABS(DY) + ABS(DX))
                        # print(X, Y)
                        pass                        
                    else : 
                        print("ICT >= 500")
                        if (IFIT == 0): 
                            if (IN < 5): 
                                (XO, YO, X, Y, ICT, IN) = POLRT_increment_initial_values(XO, YO, ICT, IN)
                            else: 
                                IER = 3
                                return (COF, ROOTR, ROOTI, IER)
                        else :
                            X = XPR
                            Y = YPR
                            IFIT = 0
                            (X, Y, N, SUMSQ, ALPHA) = POLRT_check_root(X, Y, N)
                            (COF, ROOTI, ROOTR, N2) = POLRT_compute_root(COF, ALPHA, N, SUMSQ, ROOTI, ROOTR, N2, X, Y)
                            (XO, YO, IN) = POLRT_set_initial_values(XO, YO, IN)
                            (XO, YO, X, Y, ICT, IN) = POLRT_increment_initial_values(XO, YO, ICT, IN)
                else: 
                    # (ABS(DY) + ABS(DX) < 1e-5) 
                    # print("DX+DY < 1e-5... zero found")
                    # print("COF, XCOF: ", COF, XCOF)
                    for L in range(NXX):
                        # (MT + 1) = KJ1 - (L + 1) + 1
                        MT = KJ1 - L - 1
                        # print("L, MT: ", L, MT)
                        TEMP = XCOF[MT]
                        XCOF[MT] = COF[L]
                        COF[L] = TEMP
                    # print("COF, XCOF: ", COF, XCOF)
                    ITEMP = N
                    N = NX
                    NX = ITEMP
                    if (IFIT == 0): 
                        IFIT = 1
                        XPR = X
                        YPR = Y
                    else : 
                        IFIT = 0
                        (X, Y, N, SUMSQ, ALPHA) = POLRT_check_root(X, Y, N)
                        (COF, ROOTI, ROOTR, N2) = POLRT_compute_root(COF, ALPHA, 
                            N, SUMSQ, ROOTI, ROOTR, N2, X, Y)
                        (XO, YO, IN) = POLRT_set_initial_values(XO, YO, IN)
                        (XO, YO, X, Y, ICT, IN) = POLRT_increment_initial_values(XO, YO, ICT, IN)
    # end while (N > 0)
    return (COF, ROOTR, ROOTI, IER)
#_______________________________________________________________________________
#_______________________________________________________________________________
def POLRT_set_initial_values(XO, YO, IN):
    """
    (XO, YO, IN) = POLRT_set_initial_values(XO, YO, IN)
    """
    # set initial values
    XO = .00500101
    YO = 0.01000101
    # zero initial value counter
    IN = 0
    return (XO, YO, IN)
#_______________________________________________________________________________
#_______________________________________________________________________________
def POLRT_eval_pol_and_deriv(COF, X, Y, SUMSQ, N):
    # evaluate polynomial and derivatives    
    UX = 0.0
    UY = 0.0
    V = 0.0
    YT = 0.0
    XT = 1.0
    # U = COF[N + 1]
    # print("N, COF[N]", N, COF[N])
    flag = 1
    U = COF[N]
    if (U == 0): 
        # print("U == 0")
        flag = 0
        return (flag, X, Y, SUMSQ, N, 0, 0)
    else : 
        for I in range(N): 
            # (L + 1) = N - (I + 1) + 1
            L = N - I - 1
            XT2 = X * XT - Y * YT
            YT2 = X * YT + Y * XT
            # print((X, Y, XT, YT, XT2, YT2, U, V, UX, UY))
            U += COF[L] * XT2
            V += COF[L] * YT2
            FI = (I + 1)
            UX += FI * XT * COF[L]
            UY -= FI * YT * COF[L]
            XT = XT2
            YT = YT2
        SUMSQ = UX * UX + UY * UY
        # print(SUMSQ)
        if (SUMSQ != 0): 
            DX = (V * UY - U * UX) / SUMSQ
            X += DX
            DY = -(U * UY + V * UX) / SUMSQ
            Y += DY
        else :
            DX = 0.
            DY = 0.
        # print(X, Y)
    return (flag, X, Y, SUMSQ, N, DX, DY)
#_______________________________________________________________________________
#_______________________________________________________________________________
def POLRT_increment_initial_values(XO, YO, ICT, IN): 
    """
    """
    X = XO
    # increment initial values and counter
    XO = -10.0 * YO
    YO = -10.0 * X
    # set X and Y to current value
    X = XO
    Y = YO
    ICT = 0
    IN += 1
    return (XO, YO, X, Y, ICT, IN)
#_______________________________________________________________________________
#_______________________________________________________________________________
def POLRT_check_root(X, Y, N):
    """
    if (abs(Y/X) < 1e-4) this emplies that the root is approximatley a pure real. 
    if not, we have a complex value and the solution have a complex conjugate. 
    """
    if (abs(Y / X) >= 1e-4): 
        ALPHA = X + X
        SUMSQ = X * X + Y * Y
        N -= 2
        # print("imaginary roots -> N -= 2")
    else:   
        Y = 0.0
        SUMSQ = 0.0
        ALPHA = X
        N -= 1
        # print("real root -> N -= 1")
    return (X, Y, N, SUMSQ, ALPHA)
#_______________________________________________________________________________
#_______________________________________________________________________________
def POLRT_compute_root(COF, ALPHA, N, SUMSQ, ROOTI, ROOTR, N2, X, Y): 
    # 140
    # print("compute root COF: ", COF)
    COF[1] += ALPHA * COF[0]    
    for L in range(1, N): 
        COF[L + 1] += ALPHA * COF[L] - SUMSQ * COF[L - 1]
    # print("compute root COF: ", COF)
    ROOTI[N2] = Y
    ROOTR[N2] = X
    N2 += 1
    if (SUMSQ != 0):
        Y = -Y
        SUMSQ = 0.0
        ROOTI[N2] = Y
        ROOTR[N2] = X
        N2 += 1
    # print("ROOT", ROOTR, ROOTI)
    return (COF, ROOTI, ROOTR, N2)
#_______________________________________________________________________________
#_______________________________________________________________________________
def COHERE(L, N, S): 
    """
    COHERE: COHEREnce
    
    p. 215
    
    """
    C = numpy.empty((L * N * N))
    for JP in range(1, N): 
        J = JP - 1
        for K in range(JP, N): 
            for I in range(L):
                IJK = L * N * K + L * J + I
                IKJ = L * N * J + L * K + I
                IJJ = L * N * J + L * J + I
                IKK = L * N * K + L * K + I
                num = S[IJK] ** 2. + S[IKJ] ** 2.
                den = S[IJJ] * S[IKK]
#                 print(I, J, K)
#                 print(IJK, IKJ, IJJ, IKK)
                # print(S[IJK], S[IKJ], S[IJJ], S[IKK])
                CO = (abs(num / den)) ** 0.5
                PH = ATAN2(S[IKJ], S[IJK])
#                 print(S[IKJ], S[IJK], num, den, CO, PH)
                # PH = ATAN2(S[IJK], S[IKJ])
                C[IJK] = CO
                C[IKJ] = 180.0 * (PH / pi)
    for J in range(N):
        OneJJ = L * N * J + L * J
#        print(OneJJ)
#        print(C)
        # print(S)
        # print(C)
        C = MOVE(L, S, OneJJ, C, OneJJ)
        # print(C)
#        print(C[OneJJ:])
        C[OneJJ:] = NORMAG(L, C[OneJJ:])
#        print(C)
    return C 
#_______________________________________________________________________________
#_______________________________________________________________________________
def SIMEQ1(M, N, A, B, C): 
    """
    SIMEQ1: SIMultaneous EQuation
    
    p. 39
    """
    # NMAX = 25
    # S = numpy.empty((NMAX*NMAX))
    # S = MOVE(N, N, B, S)
    # MAINE(N, S, B)
#    try : 
#        B = numpy.linalg.inv(B.reshape(N, N))
#    except numpy.linalg.linalg.LinAlgError : 
#        B = numpy.zeros((N * N,) )
#    B = B.reshape((N * N, ))
#    print(B)
    S = numpy.empty((N*N, ))
    S = MAINE(N, B, S)
#    print(S)
    A = numpy.zeros((M * N, ))
    for I in range(M): 
        for J in range(N):
            IJ = J * M + I
            A[IJ] = 0.0
            for K in range(N): 
                IK = K * M + I
                KJ = J * N + K
#                print((I, K, J, IK, KJ, IJ))
                A[IJ] += C[IK] * S[KJ]
#    B = MOVE(N * N, S, B)
#    print(A.shape)
    return A
#_______________________________________________________________________________
#_______________________________________________________________________________
def NORMEN(LX, X): 
    """
    NORMEN normalize an array by dividing each element by the RMS energy
    
    p.23
    """
    E = 0.
    for I in range(LX):
        E += X[I] * X[I]
    E = SQRT(E)
    for I in range(LX): 
        X[I] /= E
    return X
#_______________________________________________________________________________
#_______________________________________________________________________________
def NLOGN(N, X, SIGN): 
    """
    NLOGN Fast Fourier transform method 
    
    p. 63
    
    """
    NMAX = 25
    M = numpy.empty(NMAX)
    X = X.astype("complex")
    LX = 2 ** N
    for I in range(N): 
        M[I] = 2 ** (N - I - 1)
    for L in range(N): 
        NBLOCK = 2 ** L
        LBLOCK = LX // NBLOCK
        LBHALF = LBLOCK // 2
#        print("LBLOCK", LBLOCK, "LBHALF", LBHALF)
        K = 0
        for IBLOCK in range(NBLOCK): 
            FK = K 
            FLX = LX
            V = SIGN * 2 * PI * FK / FLX
            WK = complex(COS(V), SIN(V))
#            print("FK", FK, "WK", WK)
            ISTART = LBLOCK * IBLOCK
            for I in range(LBHALF): 
                J = ISTART + I
                JH = J + LBHALF
#                print("J", J, "JH", JH, "X[JH]", X[JH])
                Q = X[JH] * WK
                X[JH] = X[J] - Q
                X[J] = X[J] + Q
            for I in range(1, N): 
                II = I
                if (K < M[I]): 
                    break
                else: 
                    K -= M[I]
            K += M[II]
#       end IBLOCK
#   end L
    K = 0
    for J in range(LX):
        if (K >= J): 
            HOLD = X[J]
            X[J] = X[K]
            X[K] = HOLD
#            print("K", K, "X", X)
        for I in range(N): 
            II = I
            if (K < M[I]):
                break
            else: 
                K -= M[I]
        K += M[II]
    # end J
    if (SIGN < 0.0): 
        return X
    for I in range(LX): 
        X[I] /= FLX
    return X
#_______________________________________________________________________________
#_______________________________________________________________________________
def NORM1(LX, X): 
    """
    NORM1 normalize an array by its first element
    
    p.23
    """
    X1 = X[0]
    if (LX <= 0): 
        return X
    for I in range(LX): 
        X[I] /= X1
    return X
#_______________________________________________________________________________
#_______________________________________________________________________________
def NORMEQ(N, M, LF, R, G): 
    """
    NORMEQ: normal equation
    
    p. 247
    """
    AP = numpy.empty((N * N * LF, ))
    BP = numpy.empty((N * N * LF, ))
    VA = numpy.empty((N * N, ))
    VB = numpy.empty((N * N, ))
    DA = numpy.empty((N * N, ))
    DB = numpy.empty((N * N, ))
    CA = numpy.empty((N * N, ))
    CB = numpy.empty((N * N, ))
    CF = numpy.empty((M * N, ))
    GAM = numpy.empty((M * N, ))
    A = numpy.zeros((N * N * LF, ))
    B = numpy.zeros((N * N * LF, ))
    F = numpy.zeros((M * N * LF, ))
    for I in range(N): 
        for J in range(N):
            IJ = N * J + I
            IJ1 = N * J + I
            VA[IJ] = R[IJ1]
            VB[IJ] = R[IJ1]
        II1 = N * I + I
        A[II1] = 1.
        B[II1] = 1.
#    print(("393", F.shape, M, N))
    F[:M*N] = SIMEQ1(M, N, F, R[:N*N], G[:N*N])
#    print(("395", F.shape))
#    print("389")
    if (LF == 1): 
        return (A, B, F, AP, BP, VA, VB, DA, DB, CA, CB, CF, GAM)
#    print("392")
#    print(A, B, F, VA, VB)
    for L in range(1, LF): 
        # ZERO(N*N, DA)
        DA = numpy.zeros((N * N, ))
        OneOneL = L * M * N
        # MOVE(M * N, G[OneOneL], GAM[0])
        GAM = MOVE(M * N, G, OneOneL, GAM, 0)
        for I in range(N): 
            for LI in range(L + 1): 
                LD = L - LI
#                print("L, LI, LD", L, LI, LD)
                for K in range(N): 
                    IKLI = LI * N * N + K * N + I 
                    KILD = LD * N * N + I * N + K
                    for J in range(N): 
                        KJLD = LD * N * N + J * N + K
                        IJ = J * N + I
#                        print(KJLD, IJ)
                        DA[IJ] -= A[IKLI] * R[KJLD]
                    for J in range(M): 
                        JI = I * M + J
                        JKLI = LI * M * N + K * M + J
#                        print(JI, JKLI, KILD, GAM.shape, F.shape, R.shape)
                        GAM[JI] -= F[JKLI] * R[KILD]
            for J in range(N): 
                IJ = J * N + I
                JI = I * N + J
                DB[JI] = DA[IJ]
        CA[:N*N] = SIMEQ1(N, N, CA, VB, DA)
        CB[:N*N] = SIMEQ1(N, N, CB, VA, DB)
        AP = MOVE(N * N * (L + 1), A, 0, AP, 0)
        BP = MOVE(N * N * (L + 1), B, 0, BP, 0)
#        print(("CA, CB, AP, BP", CA, CB, AP, BP))
        for J in range(N): 
            for K in range(N): 
                KJ = J * N + K
                for LI in range(L + 1):
                    LD = L - LI
                    KJLD = LD * N * N  + J * N + K
                    for I in range(N): 
                        IJLI = LI * N * N + J * N + I
                        IK = K * N + I
                        A[IJLI] += CA[IK] * BP[KJLD]
                        B[IJLI] += CB[IK] * AP[KJLD]
                for I in range(N):
                    IJ = J * N + I
                    IK = K * N + I 
                    VA[IJ] -= CA[IK] * DB[KJ]
                    VB[IJ] -= CB[IK] * DA[KJ]
#        SIMEQ1(M, N, CF, VB, GAM)
#        print(M, N, CF, VB, GAM)
        CF = SIMEQ1(M, N, CF, VB, GAM)
        for LI in range(L + 1):
            LD = L - LI
            for J in range(N):
                for K in range(N):
                    KJLD = LD * N * N + J * N + K
                    for I in range(M): 
                        IJLI = LI * N * M + J * M + I
                        IK = K * M + I
#                        print(IJLI, IK, KJLD, F.shape)
                        F[IJLI] += CF[IK] * B[KJLD]
        # end fo LI
    # end for L
    return (A, B, F, AP, BP, VA, VB, DA, DB, CA, CB, CF, GAM)
#_______________________________________________________________________________
#_______________________________________________________________________________
def MAINV(N, A, B): 
    """
    Matrix Inversion
    
    p. 42
    """
    DET = 0
    ADJUG = numpy.zeros(N * N)
    P = numpy.zeros(N)
    (B, DET, ADJUG, P) = FADDEJ(N, A, B, DET, ADJUG, P)
    return B
#_______________________________________________________________________________
#_______________________________________________________________________________
def MAINV_CPLX(N, A, B): 
    """
    Matrix Inversion for complex number
    
    p. 42
    """
    DET = 0
    ADJUG = numpy.zeros(N * N)
    P = numpy.zeros(N)
    (B, DET, ADJUG, P) = FADDEJ_CPLX(N, A, B, DET, ADJUG, P)
    return B
#_______________________________________________________________________________
#_______________________________________________________________________________
def MAINE(N, A, B): 
    """
    MAINE: MAtrix INversion by Escalator method (symmetric matrix)
    
    p. 11
    """
    B[0] = 1.0 / A[0]
    if (N==1): 
        return B
    NN = N * N
    for I in range(1, NN):
        B[I] = 0.0
    for M in range(1, N):
        K = M 
        MM = M + M * N
        EK = A[MM]
#        print(EK)
        for I in range(K): 
            for J in range(K):
                MI = M + I * N
                IJ = I + J * N
                JM = J + M * N
                EK -= A[MI] * B[IJ] * A[JM]
#                print("EK", EK)
        B[MM] = 1.0 / EK
        for I in range(K): 
            IM = I + M * N
            for J in range(K): 
                IJ = I + J * N
                JM = J + M * N
                B[IM] -= B[IJ] * A[JM] / EK
            MI = M + I * N
            B[MI] = B[IM]
#            print(B)
        for I in range(K): 
            IM = I + M * N
            for J in range(K): 
                MJ = M + J * N
                IJ = I + J * N
                B[IJ] += B[IM] * B[MJ] * EK
#    print(M, B)
    return B 
#_______________________________________________________________________________
#_______________________________________________________________________________
def MINSN(LX, X, XMIN, INDEX): 
    """
    MINSN find the minimum element of an array, taking into account the algebraic signs of the elements. 
    
    p. 21
    """
    INDEX = 0
    for I in range(LX): 
        if (X[INDEX] > X[I]): 
            INDEX = I
        XMIN = X[INDEX]
    return (XMIN, INDEX)
#_______________________________________________________________________________
#_______________________________________________________________________________
def MAXSN(LX, X, XMAX, INDEX): 
    """
    MAXSN find the maximum element of an array, taking into account the algebraic signs of the elements. 
    
    p. 21
    """
    INDEX = 0
    for I in range(LX): 
        if (X[INDEX] < X[I]): 
            INDEX = I
        XMAX = X[INDEX]
    return (XMAX, INDEX)
#_______________________________________________________________________________
#_______________________________________________________________________________
def TRIANG(N, TOP, S, SPACE): 
    """
    TRIANG factor a positive definite symetric matrix a into the product 
    $$ a = s^t \, s$$
    where $s$ is a triangular matrix, and $s^t$ is the transpose os $s$. 
    
    p.47
    
    see also numpy.linalg.cholesky 
    """
    SPACE = ZERO(N, SPACE)
    S = ZERO(N * N, S)
    for I in range(N): 
#        IP = (I + 1) + 1 - 1
        IP = I + 1
#        IM = (I + 1) - 1, not an index just a counter
        IM = I
        II = I + I * N
        S[II] = SQRT(TOP[II] - SPACE[I])
        if (I == N - 1): 
            break
        for K in range(IP, N): 
            E = 0.
            if (IM == 0):
                pass
            else : 
                for L in range(IM):
                    LI = L + I * N
                    LK = L + K * N
                    E += S[LI] * S[LK]
            IK = I + K * N
            II = I + I * N
            X = (TOP[IK] - E) / S[II]
            SPACE[K] += X * X
            S[IK] = X
    return (S, SPACE)
#_______________________________________________________________________________
#_______________________________________________________________________________
def TRIANG_CPLX(N, TOP, S, SPACE): 
    """
    TRIANG_CPLX factor a positive definite symetric matrix a into the product 
    $$ a = s^t \, s$$
    where $s$ is a triangular matrix, and $s^t$ is the transpose os $s$. 
    
    p.47
    
    see also numpy.linalg.cholesky 
    """
    SPACE = ZERO(N, SPACE)
    S = ZERO(N * N, S)
    for I in range(N): 
        IP = I + 1
        IM = I
        II = I + I * N
        S[II] = CSQRT(TOP[II] - SPACE[I])
        if (I == N - 1): 
            break
        for K in range(IP, N): 
            E = 0.
            if (IM == 0):
                pass
            else : 
                for L in range(IM):
                    LI = L + I * N
                    LK = L + K * N
                    E += S[LI] * S[LK]
            IK = I + K * N
            II = I + I * N
            X = (TOP[IK] - E) / S[II]
            SPACE[K] += X * X
            S[IK] = X
    return (S, SPACE)
#_______________________________________________________________________________
#_______________________________________________________________________________
def TRIG(LX, X, W): 
    """
    TRIG computes one value of Fourier transform by the sum of angles TRIGonometric formula for sine and cosine. 
    
    p. 65
    """
    COSNW = 1.
    SINNW = 0.
    SINW = SIN(W)
    COSW = COS(W)
    S = 0.0
    C = 0.0
    for I in range(LX): 
        C += COSNW * X[I]
        S += SINNW * X[I]
        T = COSW * COSNW - SINW * SINNW
        SINNW = COSW * SINNW + SINW * COSNW
        COSNW = T
    return (S, C)
#_______________________________________________________________________________
#_______________________________________________________________________________
def BRAINY(NRA, NCA, LA, A, NRB, NCB, LB, B, C):
    """
    BRAINY performs matrix polynomial multiplication.
    It is the multichannel counterpart of FOLD in the single channel case. 
    Comes from the analogy with convolution in the brain. 
    
    p. 154
    """
    LC = LA + LB - 1
    # ZERO(NRA*NRB*LC, C)
    # print(NRA * NRB * LC)
    C = ZERO(NRA * NRB * LC, C)
    for I in range(LA):
        for J in range(LB): 
            K = I + J
            for M in range(NRA): 
                for N in range(NCB): 
                    MNK = M + N * NRA + K * NRA * NCB
                    for L in range(NCA):
                        MLI = M + L * NRA + I * NRA * NCA
                        LNJ = L + N * NCA + J * NCA * NCB
                        C[MNK] += A[MLI] * B[LNJ]
    return C
#_______________________________________________________________________________
#_______________________________________________________________________________
def RITE(NB, M, N, L, A): 
    """
    RITE 
    
    p. 55
    """
    i0 = 0
    while (i0 < L):
        listI = []
        nL = 0
        for i in range(NB):
            if (i0 + i < L): 
                listI.append(i0 + i)
                nL += 1
        i0 += NB
        for j in range(M):
            s = ""
            for iL in range(nL):         
                if (j == 0):
                    s += "["
                else: 
                    s += " "
                for k in range(N): 
                    if (k == 0): 
                        s += "["
                    else: 
                        s += " "
                    IJKIL = j + k * M + listI[iL] * M * N 
                    c = A[IJKIL]
                    s += str(c)
                    if k == (N - 1): 
                        s += "]"
                    else: 
                        s += " "
                if (j == M - 1):
                    s += "]"
                else: 
                    s += " "
            print(s)
        print("")
    return 

#_______________________________________________________________________________
#_______________________________________________________________________________
def SCALE(S, LX, X): 
    """
    SCALE multiply each element of an array X, of length LX, by a factor c
    
    p. 19
    """
    for I in range(LX): 
        X[I] = S * X[I]
    return X
#_______________________________________________________________________________
#_______________________________________________________________________________
def POMAIN(N, LA, A, ADJ, P, DET, S): 
    """
    POMAIN Polynomial matrix inversion
    
    p. 162
    """
    # complex A, ADJ, P, DET, S
    LADJ = (N - 1) * (LA - 1) + 1
    LDET = N * (LA - 1) + 1
    P = numpy.empty((N * LDET, ), "complex")
    DET = numpy.empty((LDET, ), "complex")
    ADJ = numpy.empty((N * N * LADJ, ), "complex")
    S = numpy.empty(N * N * LDET, "complex")
    S = MOVE(N * N * LA, A, 0, S, 0)
    J = LA
    #print("___POMAIN___")
    for L in range(N): 
        # Calculate coefficients P[., K] of characteristic polynomial
        for K in range(J):
            LK = L + K * N
            P[LK] = 0.
            for I in range(N): 
                IIK = I + I * N + K * N * N
                # P[LK] += S[IIK] / float(L)
                P[LK] += S[IIK] / float(L + 1)
        # if (L != N): 
        if (L < N - 1): 
            # Substract P[., K]*identity matrix 
            # MOVE(N * N * J, S, ADJ)
            ADJ = MOVE(N * N * J, S, 0, ADJ, 0)
            for I in range(N): 
                for K in range(J): 
                    IIK = I + I * N + K * N * N
                    LK = L + K * N
                    ADJ[IIK] -= P[LK]
            # Multiply by input matrix 
            S = BRAINY(N, N, LA, A, N, N, J, ADJ, S)
            J += LA - 1
    # Give determinant and adjugate correct sign
    # Now J = LDET = LA + (N - 1) * (LA - 1) = N * (LA - 1) + 1
    # Hence J - LA + 1 = (N - 1) * (LA - 1) + 1 = LADJ
    ADJ = SCALE(float(2 * MOD(N, 2) - 1), N * N * (J - LA + 1), ADJ)
    for L in range(J): 
        # NL = N + L * N
        NL = N - 1 + L * N
        DET[L] = P[NL] * float(2 * MOD(N, 2) - 1)
    #print("___")
    return (ADJ, P, DET, S)
#_______________________________________________________________________________
#_______________________________________________________________________________
def FACT(N, LA, A, ADJ, ZEROS, S, B): 
    """
    FACT factors the p x p matrix polynomial (1) into binomial factors (2)
    
    (1) $$ A(z) = a_0 + a_1 \, z + a_2 \, z^2+ \cdots + a_m \, z^m $$
    (2) $$ A(z) = b_0 \, (I + z \, b_1) \, (I + z \, b_2) \cdots (I + z \, b_m) $$
    
    p. 177
    
    """
    # COMPLEX A,ADJ,ZEROS,S,B,ZI
    A = A.astype("complex")
    ADJ = ADJ.astype("complex")
    ZEROS = ZEROS.astype("complex")
    S = S.astype("complex")
    B = B.astype("complex")
    ZI = complex(0, 0)
    LS = S.size / float(N * N)
    # LADJ=(N-1)*(LA-1)+1
    # DIMENSION ADJ(N,N,LADJ),ZEROS(N,LA-1),S(N,N,5+2*LADJ)
    NP = LA - 1
    # INPUTS ARE MATRIX POLYNOMIALS A,ITS ADJUGATE ADJ,ZEROS OF DET A
    # OUTPUTS IS B WHERE A(Z)=B0*(I+ZB1)*(I+ZB2)*...*(1+ZBNP)
    # MOVE ADJUGATE TO S(1, 1, 6)
    OOS = 5 * N * N 
    # S = MOVE(N * N * ((N - 1) * NP + 1), ADJ, 0, S, OOS)
    # print("ADJ.shape, S.shape, N * N * ((N - 1) * NP + 1), OOS", ADJ.shape, S.shape, N * N * ((N - 1) * NP + 1), OOS)
    S = MOVE(N * N * ((N - 1) * NP + 1), ADJ, 0, S, OOS)
    # LOOP ON THE FACTORS (I + Z * Q * D * QINV)
    for IP in range(NP):
        print("-----------------")
        print("IP ", IP)
        print(" ")
        # L = (N - 1) * NP + IP
        L = (N - 1) * NP + (IP + 1)
        # CALL RITE(2, N, N, L, S(1, 1, 6))
        S = ZERO(N * N * 5, S) 
        print("___________")
        print("S")
        RITE(4, N, N, LS, S)
        print("___________")

        for ICOL in range(N): 
            print(".......")
            print("ICOL, IP", ICOL, IP)
            # INSERT ROOTS INTO (PREVIOUS FACTORS*ADJUGATE) TO GET EIGENCOLUMNS.
            B = ZERO(N * N, B)
            for I in range(L): 
                # ZI = ZEROS[ICOL, IP] ** (I - 1)
                ICOLIP = ICOL + IP * N
                ZI = ZEROS[ICOLIP] ** I
                print("I, ZI, ", I, ZI)
                for K in range(N): 
                    for J in range(N):
                        # JKO = (J, K, 1)
                        JKO = J + K * N
                        # JKIPF = (J, K, I + 5)
                        JKIPF = JKO + (I + 5) * N * N
                        print("K,J,JKO, JKIPF, S[JKIPF], B[JKO], ZI", K,J,JKO, JKIPF, S[JKIPF].real, B[JKO], ZI)
                        B[JKO] += S[JKIPF] * ZI
            # 10
            print("___________")
            print("B")
            RITE(3, N, N, LA, B)
            print("___________")
            # print("B", B)
            # FORM MATRIX Q OF EIGENCOLUMNS AND STORE IN S(1,1,3)
            # OICOLT = (1, ICOL, 3)
            OICOLT = ICOL * N + 2 * N * N
            S = MOVE(N, B, 0, S, OICOLT)
            # FORM DIAGONAL MATRIX D WITH -1./ZERO AND STORE IN S(1,1,4)
            # ICOLICOLF = (ICOL, ICOL, 4)
            ICOLICOLF = ICOL + ICOL * N + 3 * N * N
            # ICOLIP = (ICOL, IP)
            ICOLIP = ICOL + IP * N
            print("ICOLICOLF", ICOLICOLF)
            S[ICOLICOLF] = -1. / ZEROS[ICOLIP]
            # FORM IDENTITY MATRIX I AND STORE IN S(1,1,1).
            # ICOLICOLO = (ICOL, ICOL, 1) 
            ICOLICOLO = ICOL + ICOL * N 
            print("ICOLICOLO", ICOLICOLO)
            S[ICOLICOLO] = 1.0
            print("___________")
            print("S")
            RITE(4, N, N, LS, S)
            print("___________")
            print(".......")
        # 20 
        # FORM Q**-1 AND STORE IN S(1,1,5).
        # OOTW = (1, 1, 2)
        OOTW = 1 * N * N
        # OOTH = (1, 1, 3)
        OOTH = 2 * N * N
        # OOFO = (1, 1, 4)
        OOFO = 3 * N * N
        # OOFI = (1, 1, 5)
        OOFI = 4 * N * N
        S[OOFI: ] = MAINV_CPLX(N, S[OOTH: ], S[OOFI: ])
        print("______________")
        print("S")
        RITE(4, N, N, LS, S)
        print("______________")
       
        # FORM Q*D*Q**-1 AND STORE IN S(1,1,2)
        print("______________")
        print("Q")
        RITE(1, N, N, 1, S[OOTH: ])
        print("______________")
        print("______________")
        print("D")
        RITE(1, N, N, 1, S[OOFO: ])
        print("______________")
        B[: N * N] = BRAINY(N, N, 1, S[OOTH: ], N, N, 1, S[OOFO: ], B[: N * N])
        print("______________")
        print("B apres Q D")
        RITE(3, N, N, 3, B)
        print("______________")
        print("______________")
        print("Q-1")
        RITE(1, N, N, 1, S[OOFI: ])
        print("______________")
        
        S[OOTW: OOTW + N * N] = BRAINY(N, N, 1, B, N, N, 1, S[OOFI: ], S[OOTW: OOTW + N * N])
        print("______________")
        print("S apres Q D Q-1")
        RITE(4, N, N, LS, S)
        print("______________")
        # print("S[OOTW: ]", S[OOTW: ])
        # STORE MATRIX Q*D*Q**-1 IN B ARRAY.
        # I = NP - IP + 2
        I = NP - IP
        OOI = I * N * N
        print("I, OOTW, OOI", I, OOTW, OOI)
        B = MOVE(N * N, S, OOTW, B, OOI)
        print("______________")
        print("B")
        RITE(3, N, N, 3, B)
        print("______________")
        # NEW FACTOR = (I+Z*Q*D*Q**-1) = (I+Z*B)
        # FORM FACTOR*(PREVIOUS FACTORS*ADJUGATE) AND STORE IN S(1,1,6)
        # OOLPS = (1, 1, L + 6)
        OOLPS = (L + 5) * N * N
        # OOS = (1, 1, 6)
        OOS = 5 * N * N
        print("L, S.shape, OOLPS, OOS", L, S.shape, OOLPS, OOS)
        # print(OOLPS + N * N * (L + 1))
        S[OOLPS: OOLPS + N * N * (L + 1)] = BRAINY(N, N, 2, S, N, N, L, S[OOS: ], S[OOLPS: OOLPS + N * N * (L + 1)])
        print("____")
        print("S avant MOVE")
        RITE(4, N, N, LS, S)
        print("____")
        S = MOVE(N * N * (L + 1), S, OOLPS, S, OOS)
        print("____")
        print("S")
        RITE(4, N, N, LS, S)
        print("____")
        print("-------------")
    # 30 
    B = MOVE(N * N, A, 0, B, 0)
    return (S, B)
#_______________________________________________________________________________
#_______________________________________________________________________________
def FACTOR(M, N, BS): 
    """
    FACTOR
    
    p. 197
    
    """
    NX = 49
    XCOF = numpy.zeros((NX + 1, ))
    COF = numpy.zeros((NX + 1, ))
    RZR = numpy.zeros((NX, ))    
    CZR = numpy.zeros((NX, ))
    N1 = N - 1
    NZB = N1 * M
    LR = N + N1
    LAJ = (M - 1) * (LR - 1) + 1
    LDETR = (LR - 1) * M + 1
    NZR = 2 * NZB
    LRR = LAJ + N1
    print("LR, LAJ, LDETR, LRR : ", LR, LAJ, LDETR, LRR)
    # initialement R = numpy.zeros((M * M * LR, ), "complex")
    R = numpy.zeros((M * M * LR, ), "complex")
    AJ = numpy.zeros((M * M * LAJ), "complex")
    DETR = numpy.zeros((LDETR, ), "complex")
    DET = complex(0, 0)
    DETP = complex(0, 0)
    ZR = numpy.zeros((NZR, ), "complex")
    ZB = numpy.zeros((M * (N - 1), ), "complex")
    RR = numpy.zeros((M * M * LRR, ), "complex")
    FACT = numpy.zeros((M * M * 2, ), "complex")
    DIAG = numpy.zeros((M * M, ), "complex")
    CR = numpy.zeros((M * M, ), "complex")
    P = numpy.zeros((M * M), "complex")
    PINV = numpy.zeros((M * M, ), "complex")
    V = numpy.zeros((M * M, ), "complex")
    TEMP = numpy.zeros((M * M * LDETR, ), "complex")
    B = numpy.zeros((M * M * N, ), "complex")
    BZ1 = numpy.zeros((M * M, ), "complex")
    RZ1 = numpy.zeros((M * M, ), "complex")
    W = numpy.zeros((M * M, ), "complex")
    BZINV = numpy.zeros((M * M, ), "complex")
    # complex DETP, DET, ZER, Z1
    # Compute one side of autocorrelation
    OON = 0 + 0 * M + (N - 1) * M * M  
    R[OON: OON + N * M * M] = HEAT_CPLX(M, M, N, BS, M, M, N, BS, N, R[OON: ])
    # Generate other side of autocorrelation
    # print(("R", R))
    for I in range(1, N): 
        # (JP + 1) = N + (I + 1) - 1
        JP = N + I - 1 
        # (JM + 1) = N - (I + 1) + 1
        JM = N - I - 1
        for J in range(M): 
            for K in range(M): 
                JKJM = J + K * M + JM * M * M
                KJJP = K + J * M + JP * M * M
                # print("I, JP, JM, J, K, JKJM, KJJP", I, JP, JM, J, K, JKJM, KJJP)
                R[JKJM] = R[KJJP]
    # print("R", R)
    # Find adjugate and autocorrelation of R
    # print("AJ.shape, RR.shape, DETR.shape, TEMP.shape", AJ.shape, RR.shape, DETR.shape, TEMP.shape)
    LPOMAINRR = M * (M * (LR - 1) + 1)
    (AJ, RR[: LPOMAINRR], DETR, TEMP) = POMAIN(M, LR, R, AJ, RR, DETR, TEMP)
    # print("AJ.shape, RR.shape, DETR.shape, TEMP.shape", AJ.shape, RR.shape, DETR.shape, TEMP.shape)
    # Find zeros of determinant
    for I in range(LDETR): 
        XCOF[I] = DETR[I].real
    IER = 0
    (COF, RZR, CZR, IER) = POLRT(XCOF, COF, LDETR - 1, RZR, CZR, IER)
    # print(RZR, CZR, XCOF, LDETR - 1)
    for I in range(NZR): 
        ZR[I] = complex(RZR[I], CZR[I])
    J = 0
    # print(ZR)
    for I in range(NZR): 
        # It is assumed that no zero has magnitude unity
        if (numpy.abs(ZR[I]) >= 1.0): 
            # CABS is abs for complex 
            # Choose the zeros with magnitude greater than unity to form B
            DETR[J] = ZR[I]
            J += 1
    # print(ZR, DETR)
    ZB = MOVE(NZB, DETR, 0, ZB, 0)
    # Set B = I, FACT(J, K, 1) = 1        
    FACT = ZERO(M * M, FACT)
    B = ZERO(M * M, B)
    for J in range(M):
        JJ0 = J + J * M 
        FACT[JJ0] = 1. 
        B[JJ0] = 1. 
    LBT = 1
    # Set RR = AJ
    RR = MOVE(LAJ * M * M, AJ, 0, RR, 0)
    LRRT = LAJ
    # Loop on the factors (I - Z * PINV * DIAG * P)
    for IFACTS in range(N1): 
        # print("IFACTS, N1", IFACTS, N1)
        # Form diagonal matrix
        # ZERO(M * M, DIAG)
        DIAG = ZERO(M * M, DIAG)
        for I in range(M):
            II = I + I * M
            IIFACTS = I + IFACTS * M
            # print("IIFACTS, I, IFACTS:", IIFACTS, I, IFACTS)
            DIAG[II] = 1. / ZB[IIFACTS]
        # Insert zeros in RR to get eigenvectors
        for IVECT in range(M): 
            # ZERO(M * M, CR)
            CR = ZERO(M * M, CR)
            IVECTIFACTS = IVECT + IFACTS * M
            # print("IVECTIFACTS, IVECT, IFACTS:", IVECTIFACTS, IVECT, IFACTS)
            ZER = ZB[IVECTIFACTS]
            for I in range(LRRT):
                # ZI = ZER ** ((I + 1) - 1 - LAJ / 2)
                ZI = ZER ** (I - LAJ / 2)
                for K in range(M):
                    for J in range(M): 
                        JK = J + K * M
                        JKI = J + K * M + I * M * M
                        # print(I, J, K, JK, JKI, RR.shape, CR.shape)
                        CR[JK] += RR[JKI] * ZI 
            for I in range(M):
                IVECTI = IVECT + I * M
                # [0, I] becomes [I * M]
                P[IVECTI] = CR[I * M]
        # end for IVECT
        # Form PINV * DIAG * P
        # [1, 1, 2] becomes [1*M*M]
        OOT = 0 + 0 * M + 1 * M * M
        # print(OOT, FACT.real, DIAG)
        FACT[OOT: OOT + M * M] = BRAINY(M, M, 1, DIAG, M, M, 1, P, FACT[OOT:])
        # print(PINV.shape, DETP, DIAG.shape, TEMP.shape)
        # print(PINV)
        # try : 
        (PINV, DETP, DIAG, TEMP) = FADDEJ_CPLX(M, P, PINV, DETP, DIAG, TEMP)
        # except ComplexWarning as e:
        # print(e)
        # print(PINV)
        # print(PINV.shape, DETP, DIAG.shape, TEMP.shape)
        P = MOVE(M * M, FACT, OOT, P, 0)
        FACT[OOT: OOT + M * M] = BRAINY(M, M, 1, PINV, M, M, 1, P, FACT[OOT:])
        FACT[OOT: OOT + M * M] = SCALE(-1., M * M, FACT[OOT:])
        # Set RR = RR * (I - Z * PINV * DIAG * P)
        # print(TEMP.shape)
        # ltemp = M * M * (LRRT + 2 - 1)
        TEMP = BRAINY(M, M, LRRT, RR, M, M, 2, FACT, TEMP)
        LRRT += 1
        RR = MOVE(LRRT * M * M, TEMP, 0, RR, 0)
        # Set B = B * (I - FACT)
        # ltemp = M * M * (LBT + 2 - 1)
        TEMP = BRAINY(M, M, LBT, B, M, M, 2, FACT, TEMP)
        LBT += 1
        B = MOVE(M * M * LBT, TEMP, 0,  B, 0)
    # 60 
    TEMP = HEAT_CPLX(M, M, N, B, M, M, N, B, N, TEMP)
    # Form B(Z = 1) and R(Z = 1)
    RZ1 = ZERO(M * M, RZ1)
    BZ1 = ZERO(M * M, BZ1)
    for J in range(M): 
        for K in range(M): 
            JK = J + K * M
            for I in range(N):
                JKI = JK + I * M * M
                BZ1[JK] += B[JKI]
            for I in range(LR): 
                JKI = JK + I * M * M
                RZ1[JK] += R[JKI]
    # end for J, K, I
    # Form BINV(Z = 1) * R(Z = 1) * BINVTRANSP(Z = 1) = W
    (BZINV, DET, TEMP, W) = FADDEJ_CPLX(M, BZ1, BZINV, DET, TEMP, W)
    TEMP = BRAINY(M, M, 1, BZINV, M, M, 1, RZ1, TEMP)
    W = HEAT_CPLX(M, M, 1, TEMP, M, M, 1, BZINV, 1, W)
    # Triangularize W
    #print("V, W, TEMP", V, W, TEMP)
    (V, TEMP) = TRIANG_CPLX(M, W, V, TEMP)
    #print("V, W, TEMP", V, W, TEMP)
    # Put B in causal-chain form
    TEMP = HEAT_CPLX(M, M, N, B, M, M, 1, V, N, TEMP)
    B = MOVE(M * M * N, TEMP, 0, B, 0)
    TEMP = HEAT_CPLX(M, M, N, B, M, M, N, B, N, TEMP)
    # return (LR, R, LAJ, AJ, LDETR, DETR, NZR, ZR, ZB, LRR, RR, FACT, DIAGCR, P, PINV, V, TEMP, B, BZ1, RZ1, W, BZINV)
    return (LR, R, LAJ, AJ, LDETR, DETR, NZR, ZR, ZB, B)
#_______________________________________________________________________________
#_______________________________________________________________________________
def FADDEJ(N, A, AINV, DET, ADJUG, P): 
    """
    FADDEJ inverts a (not necessarily symmetric) n x n matrix a by a method given by Faddeev and Sominskii. 
    
    p. 40
    """
    AINV = MOVE(N * N, A, 0, AINV, 0)
    for K in range(N): 
        P[K] = 0.0
        for I in range(N): 
            II = I + I * N
            P[K] += AINV[II]
        P[K] /= float(K + 1)
        if (K == N - 1): 
            break
        ADJUG = MOVE(N * N, AINV, 0, ADJUG, 0)
        for I in range(N): 
            II = I + I * N
            ADJUG[II] = AINV[II] - P[K]
        AINV = BRAINY(N, N, 1, A, N, N, 1, ADJUG, AINV)
    AINV = MOVE(N * N, ADJUG, 0, AINV, 0)
    if (ABS(P[N - 1]) >= 1e-30): 
        AINV = SCALE(1. / P[N - 1], N * N, AINV)
    DET = P[N - 1]
    if (MOD(N, 2) == 1): 
        pass
    else : 
        DET = -DET
        for I in range(N): 
            for J in range(N): 
                IJ = I + J * N
                ADJUG[IJ] = - ADJUG[IJ]
    return (AINV, DET, ADJUG, P)
#_______________________________________________________________________________
#_______________________________________________________________________________
def FADDEJ_CPLX(N, A, AINV, DET, ADJUG, P): 
    """
    FADDEJ inverts a (not necessarily symmetric) n x n matrix a by a method given by Faddeev and Sominskii. 
    
    p. 40
    """
    # print("__FADDEJ_CPLX__")
    # print(A.dtype, AINV.dtype, DET, ADJUG.dtype, P.dtype)
    A = A.astype("complex")
    AINV = AINV.astype("complex")
    ADJUG = ADJUG.astype("complex")
    P = P.astype("complex")
    DET = complex(0, 0)
    AINV = MOVE(N * N, A, 0, AINV, 0)
    for K in range(N): 
        P[K] = 0.0
        for I in range(N): 
            II = I + I * N
            P[K] += AINV[II]
        P[K] /= float(K + 1)
        if (K == N - 1): 
            break
        ADJUG = MOVE(N * N, AINV, 0, ADJUG, 0)
        for I in range(N): 
            II = I + I * N
            ADJUG[II] = AINV[II] - P[K]
        AINV = BRAINY(N, N, 1, A, N, N, 1, ADJUG, AINV)
    AINV = MOVE(N * N, ADJUG, 0, AINV, 0)
    if (CABS(P[N - 1]) >= 1e-30):
        for I in range(N):
            for J in range(N): 
                IJ = I + J * N
                AINV[IJ] = AINV[IJ] / P[N - 1]
    DET = P[N - 1]
    if (MOD(N, 2) == 1): 
        pass
    else : 
        DET = -DET
        for I in range(N): 
            for J in range(N): 
                IJ = I + J * N
                ADJUG[IJ] = - ADJUG[IJ]
    # print("___")
    return (AINV, DET, ADJUG, P)
#_______________________________________________________________________________
#_______________________________________________________________________________
def FOLD(LA, A, LB, B, LC, C):
    """
    FOLD performs polynomial multiplication or equivalently the complete transient convolution of two signals. 
    
    p.29
    """
    LC = LA + LB - 1
    C = numpy.zeros((LC, ))
    for I in range(LA): 
        for J in range(LB): 
            K = I + J
            C[K] += A[I] * B[J] 
    return (LC, C)
#_______________________________________________________________________________
#_______________________________________________________________________________
def plotMacro(L, N, G): 
    for i in range(N):
        for j in range(N): 
            k = N * i + j + 1
            OneIJ = L * N * j + L * i 
            plt.subplot(N, N, k)
            plt.plot(G[OneIJ: OneIJ+L])
            plt.title("$\phi$(" + str(i) + ", " + str(j) + ")") 
    plt.tight_layout()
    return    
#_______________________________________________________________________________
#_______________________________________________________________________________
def plotCohere(L, N, C, Fs=0):
    if (Fs == 0): 
        f = numpy.arange(0, L)
    else : 
        f = numpy.arange(0, Fs / 2. + Fs / (2. * (L - 1)), Fs / (2. * (L - 1)))
    for i in range(N):
        for j in range(N): 
            k = N * i + j + 1
            OneIJ = L * N * j + L * i 
            plt.subplot(N, N, k)
            if (i < j): 
                plt.plot(f, C[OneIJ: OneIJ+L])
                plt.axis(ymin=0., ymax=1.)
                plt.title("|K(" + str(i) + ", " + str(j) + ")|") 
            if (i == j): 
                plt.plot(f, C[OneIJ: OneIJ+L])
                plt.axis(ymin=0., ymax=1.)
                plt.title("$\phi$(" + str(i) + ", " + str(j) + ")") 
            if (i > j): 
                plt.plot(f, C[OneIJ: OneIJ+L])
                plt.axis(ymin=-190., ymax=190.)
                plt.title(r"$\theta$(" + str(j) + ", " + str(i) + ")") 
    plt.tight_layout()
    return
#_______________________________________________________________________________
#_______________________________________________________________________________
def unitfreq2freq(L, Fs):
    """
    L = m + 1
    """
    m = L - 1
    dFU = 1 / (2 * m)
    dF = dFU * Fs
    return dF
#_______________________________________________________________________________


