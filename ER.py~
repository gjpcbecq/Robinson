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
REAL = numpy.real
AIMAG = numpy.imag
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
def eval_pol_and_deriv(COF, X, NX, NXX, Y, SUMSQ, ALPHA, N):
    # evaluate polynomial and derivatives    
    UX = 0.0
    UY = 0.0
    V = 0.0
    YT = 0.0
    XT = 1.0
    U = COF[N]
    if (U == 0): 
        X = 0.0
        NX -= 1
        NXX -= 1
        Y = 0.0
        SUMSQ = 0.0
        ALPHA = X
        N -= 1
        return (U, X, NX, NXX, Y, SUMSQ, ALPHA, N, 0, 0)
    for I in range(N): 
        L = N - I - 1
        XT2 = X * XT - Y * YT
        YT2 = X * YT + Y * XT
        # print((X, Y, XT, YT, XT2, YT2, U, V, UX, UY))
        U += COF[L] * XT2
        V += COF[L] * YT2
        FI = I
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
    return (U, X, NX, NXX, Y, SUMSQ, ALPHA, N, DX, DY)
    
def increment_initial_values(XO, YO, ICT, IN): 
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
    
def compute_root(COF, ALPHA, N, SUMSQ, ROOTI, ROOTR, N2, X, Y): 
    # 140
    COF[1] += ALPHA * COF[0]
    for L in range(1, N): 
        COF[L + 1] += ALPHA * COF[L] - SUMSQ * COF[L - 1]
    print(N2, X, Y)
    ROOTI[N2] = Y
    ROOTR[N2] = X
    N2 += 1
    if (SUMSQ != 0):
        Y = -Y
        SUMSQ = 0.0
        ROOTI[N2] = Y
        ROOTR[N2] = X
        N2 += 1
    return (COF, ROOTI, ROOTR, N2)
    
#_______________________________________________________________________________
#_______________________________________________________________________________
def AUGURY(N, LX, X, LR, R, SPIKE, FLOOR, LF, F, LY, Y, ERROR):
    """
    AUGURY
    
    p. 99
    """
    NMAX = 9
    VF = numpy.empty((NMAX * NMAX,))
    VB = numpy.empty((NMAX * NMAX, ))
    R = HEAT(N, 1, LX, X, N, 1, LX, X, LR, R)
    RT = 0.0
    for I in range(N): 
        J = I * N + I
        print(I, J)
        R[J] *= (1. + SPIKE)
        RT += R[J]
    NNLR = N * N * LR
    for L in range(LR): 
        (F, Y, VF, VB, Y) = OMEN(N, L, R, F, Y, 
            VF, VB, Y[NNLR:])
        Q = 0.0
        for I in range(N): 
            J = I * N + I
            Q += VF[J]
        ERROR[L] = Q / RT
        LF = L
        if (ERROR[L] <= FLOOR): 
            break
    LY = LX + LF
    Y = BRAINY(N, N, LR, F, N, 1, LX, X, Y)
    return (R, LF, F, LY, Y, ERROR)
#_______________________________________________________________________________
#_______________________________________________________________________________
def OMEN(N, L, R, AF, AB, VF, VB, SP): 
    """
    OMEN 
    
    p. 99
    """
    print("AF: ", AF)
    NMAX = 9
    DF = numpy.empty((NMAX * NMAX,))
    DB = numpy.empty((NMAX * NMAX, ))
    CF = numpy.empty((NMAX * NMAX,))
    CB = numpy.empty((NMAX * NMAX, ))
    L += 1
    if (L == 1):
        VF = MOVE(N * N, R, 0, VF, 0)
        VB = MOVE(N * N, R, 0, VB, 0)
        AF[:N * N] = ZERO(N * N)
        for I in range(N):
            IIZ = I + I * N 
            AF[IIZ] = 1.
        AB = MOVE(N * N, AF, 0, AB, 0)
        print("AF, AB, VF, VB, SP", AF, AB, VF, VB, SP)
        return (AF, AB, VF, VB, SP)
    ZZO = 1
    DB = HEAT(N, N, L - 1, AB, N, N, L - 1, R[ZZO:], 1, DB)
    print(R, DB)
    for I in range(N): 
        for J in range(N): 
            IJ = I * N + J
            JI = J * N + I
            DF[IJ] = DB[JI]
    SP = MAINE(N, VB, SP)
    CF = BRAINY(N, N, 1, DF, N, N, 1, SP, CF)
    SP = MAINE(N, VF, SP)
    CB = BRAINY(N, N, 1, DB, N, N, 1, SP, CB)
    SP = MOVE(N * N * (L - 1), AB, 0, SP, ZZO)
    SP[: N * N] = ZERO(N * N)
    AB = MOVE(N * N * L, SP, 0, AB, 0)
    ZZL = L - 1
    AF[ZZL: ZZL + N * N] = ZERO(N * N)
    print("----")
    print(AF, SP, DB, DF)
    AF = FORM(N, N, L, AB, CB, AF)
    SP = FORM(N, N, L, AF, CF, SP)
    DB = FORM(N, N, 1, VF, CF, DB)
    DF = FORM(N, N, 1, VB, CB, DF)
    print(AF, SP, DB, DF)
    print("----")
    return (AF, AB, VF, VB, SP)
#_______________________________________________________________________________
#_______________________________________________________________________________
def FORM(M, N, L, A, B, C): 
    """
    FORM 
    
    p. 100
    """
    for I in range(M): 
        for J in range(N): 
            for II in range(N): 
                for K in range(L):
                    IJK = I + J * M + II * M * N
                    III = I + II * M 
                    IIJK = II + J * N + K * N * N
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
#_______________________________________________________________________________
def ZERO(LX):
    """
    if (LX <= 0): 
        return
    for I in range(LX): 
        X[I] = 0
    return
    """
    X = numpy.zeros((LX, ))
    return X
#_______________________________________________________________________________
#_______________________________________________________________________________
def MOVE(LX, X, IX, Y, IY):
    """
    MOVE: MOVE one array from one storage location to another
    
    p. 18
    
    """
    cond = (XLOCF(X) - XLOCF(Y))
    # print(cond)
    if (cond < 0):
        K = LX - 1
        for I in range(0, LX): 
            Y[IY + K] = X[IX + K]
            K -= 1
    elif (cond > 0): 
        for I in range(0, LX): 
            Y[IY + I] = X[IX + I]
    else : # (cond == 0)
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
def CROSS(LX, X, LY, Y, LG):
    """
    CROSS: CROSS correlation
    
    p. 27
    
    """
    G = numpy.zeros((LG, ))
    for J in range(LG):
        # print(MIN0((LY, LX - J)))
        # N - J samples
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
def MACRO(N, LX, X, LY, Y, LG): 
    """
    MACRO: empirical MultichAnnel CROss correlation 

    p. 203
    """
    G = numpy.zeros((LG * N * N))
    for I in range(N): 
        I1 = I * LX
        for J in range(N): 
            J1 = J * LY
            IJ = LG * I + LG * N * J 
            G[IJ:IJ+LG] = CROSS(LX, X[I1:], LY, Y[J1:], LG)
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
        A[L] = - D / V
        if (L != 1):
            L1 = (L - 1) / 2
            L2 = L1 + 1
            if (L2 >= 2): 
                for J in range(1, L2): 
                    HOLD = A[J]
                    K = L - J - 1
                    A[J] += A[L] * A[K]
                    A[K] += A[L] * HOLD
            if (2 * L1 != L - 4): 
                A[L2 + 1] += A[L] * A[L2 + 1]
        V += A[L] * D
        F[L] = (G[L] - Q) / V
        L3 = L
        for J in range(L3): 
            K = L - J
            F[J] += F[L] * A[K]
            print(J, L, K)
        print((LR, L, numpy.round(F, 2)))
        if (L == LR - 1): 
            return (F, A)
        D = 0.0
        Q = 0.0
        for I in range(L + 1):
            K = L - I + 1
            print(I, K)
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
        print("K: {0}".format(K))
        SPACE = IMPULS(LR, SPACE, K)
        print(numpy.round(SPACE, 2))
        J = LR * K
        (Q[J: J + LR], SPACE[LR: 2 * LR]) = EUREKA(LR, R, SPACE, 
            Q[J : J + LR], SPACE[LR: 2 * LR])
        print(numpy.round(Q, 2))
    return (Q, SPACE)
#_______________________________________________________________________________
#_______________________________________________________________________________
def SHAPE(LB, B, LD, D, LA, A, LC, C, ASE, SPACE):
    """
    SHAPE computes the waveshaping filter
    
    p. 75
    """
    SPACE = CROSS(LB, B, LB, B, LA, SPACE)
    SPACE[LA : LA + LA + 1] = CROSS(LD, D, LB, B, LA, SPACE)
    SPACE[2 * LA + 1: ] = EUREKA(LA, SPACE, SPACE[LA:LA + LA + 1], A, 
        SPACE[2 * LA + 1: ])
    DD = DOT(LD, D, D, DD)
    AG = DOT(LA, A, SPACE[LA: LA + 1], AG)
    ASE = (DD - AG) / DD
    C = FOLD(LA, A, LB, B, LC, C)
    return (A, LC, C, ASE, SPACE)
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
    
    p.207
    """
#    ZERO(NRX * NRY * LG, G)
    print(NRX, NRY, LG)
    G = numpy.zeros((NRX * NRY * LG, ))
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
                        G[MNJ] += X[MLK] * Y[NLI]
    return G
#_______________________________________________________________________________
#_______________________________________________________________________________
def CPLX_HEAT(NRX, NCX, LX, X, NRY, NCY, LY, Y, LG, G): 
    """
    HEAT multichannel autocorrelation or cross correlation
    
    p.207
    """
#    ZERO(NRX * NRY * LG, G)
    G = numpy.zeros((NRX * NRY * LG, ))
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
    
    not working yet
    """
    flagDoubleRoot = 0 
    ROOTR = numpy.zeros((M, ))
    ROOTI = numpy.zeros((M, ))
    IFIT = 0
    N = M 
    IER = 0
    if (XCOF[N] == 0):
# set error code to 4
        IER = 4
        return (ROOTR, ROOTI, IER)
    if (N <= 0): 
# set error code to 1
        IER = 1
        return (ROOTR, ROOTI, IER)
    if (N > 36): 
# set error code to 2
        IER = 2
        return (ROOTR, ROOTI, IER)
    NX = N
    NXX = N + 1
    N2 = 0
    KJ1 = N + 1
    ICT = 0
    SUMSQ = 0.
    ALPHA = 0.
    for L in range(KJ1): 
        MT = KJ1 - L - 1
        COF[MT] = XCOF[L]
    while (N > 0): 
        print((N, COF))
    # set initial values
        XO = .00500101
        YO = 0.01000101
    # zero initial value counter
        IN = 0
        (XO, YO, X, Y, ICT, IN) = increment_initial_values(XO, YO, ICT, IN)
        (U, X, NX, NXX, Y, SUMSQ, ALPHA, N, DX, DY) = eval_pol_and_deriv(
            COF, X, NX, NXX, Y, SUMSQ, ALPHA, N)
        while (ABS(DY) + ABS(DX) > 1e-5): 
            # step iteration counter
            ICT += 1
            # print(ABS(DY) + ABS(DX))
            if (ICT < 500): 
                (U, X, NX, NXX, Y, SUMSQ, ALPHA, N, DX, DY) = eval_pol_and_deriv(COF, X, NX, NXX, Y, SUMSQ, ALPHA, N)
            else : 
                if (IN < 5): 
                    (XO, YO, X, Y, ICT, IN) = increment_initial_values(XO, YO, ICT, IN)
                    (U, X, NX, NXX, Y, SUMSQ, ALPHA, N, DX, DY) = eval_pol_and_deriv(COF, X, NX, NXX, Y, SUMSQ, ALPHA, N)
                else: 
                    IER = 3
                    return (ROOTR, ROOTI, IER)
        if (U != 0): 
            for L in range(NXX):
                MT = KJ1 - L - 1
                TEMP = XCOF[MT]
                XCOF[MT] = COF[L]
                COF[L] = TEMP
            ITEMP = N
            N = NX
            NX = ITEMP
            if (IFIT == 0): 
                IFIT = 1
                XPR = X
                YPR = Y
                (U, X, NX, NXX, Y, SUMSQ, ALPHA, N, DX, DY) = eval_pol_and_deriv(COF, X, NX, NXX, Y, SUMSQ, ALPHA, N)
                flagConv = 1
                while (ABS(DY) + ABS(DX) > 1e-5): 
                    # step iteration counter
                    ICT += 1
                    if (ICT < 500): 
                        (U, X, NX, NXX, Y, SUMSQ, ALPHA, N, DX, DY) = eval_pol_and_deriv(COF, X, NX, NXX, Y, SUMSQ, ALPHA, N)
                    else: 
                        DX = 0
                        DY = 0
                        flagConv = 0
                if (flagConv == 1): 
                    for L in range(NXX):
                        MT = KJ1 - L - 1
                        TEMP = XCOF[MT]
                        XCOF[MT] = COF[L]
                        COF[L] = TEMP
                    print(ITEMP, N, NX)
                    ITEMP = N
                    N = NX
                    NX = ITEMP
                else : 
                    X = XPR
                    Y = YPR
                    flagConv = 0
            IFIT = 0
            if (abs(Y / X) >= 1e-4): 
                ALPHA = X + X
                SUMSQ = X * X + Y * Y
                N -= 2
                flagDoubleRoot = 1
            else: #(flagDoubleRoot == 0):  
                Y = 0.0
                SUMSQ = 0.0
                ALPHA = X
                N -= 1
                flagDoubleRoot = 0
        (COF, ROOTI, ROOTR, N2) = compute_root(COF, ALPHA, N, SUMSQ, ROOTI, ROOTR, N2, X, Y)
        print(N, COF)
    # end while (N > 0)
    return 
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
#                print(0, 0, EK)
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
        if (X[INDEX] < X[I]): 
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
        if (X[INDEX] > X[I]): 
            INDEX = I
        XMAX = X[INDEX]
    return (XMAX, INDEX)
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
# ZERO(NRA*NRB*LC,C)
    C = numpy.zeros((NRA * NRB * LC, ))
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
#    complex A, ADJ, P, DET, S
    LADJ = (N - 1) * (LA - 1) + 1
    LDET = N * (LA - 1) + 1
    P = numpy.empty((N * LDET, ), "complex")
    DET = numpy.empty((LDET, ), "complex")
    ADJ = numpy.empty((N * N * LADJ, ), "complex")
#    MOVE(N * N * LA, A, S)
    S = A.copy()
    J = LA
    for L in range(N): 
# Calculate coefficients P[., K] of characteristic polynomial
        for K in range(J):
            LK = L + K * N
            P[LK] = 0.
            for I in range(N): 
                IIK = I + I * N + K * N * N
#                P[LK] += S[IIK] / float(L)
                P[LK] += S[IIK] / float(L + 1)
#        if (L != N): 
        if (L < N - 1): 
# Substract P[., K]*identity matrix 
#            MOVE(N * N * J, S, ADJ)
            ADJ = MOVE(N * N * J, S, 0, ADJ, 0)
            for I in range(N): 
                for K in range(J): 
                    IIK = I + I * N + K * N * N
                    LK = L + K * N
                    ADJ[IIK] -= P[LK]
# Multiply by input matrix 
#            BRAINY(N, N, LA, A, N, N, J, ADJ, S)
            S = BRAINY(N, N, LA, A, N, N, J, ADJ, S)
            J += LA - 1
# Give determinant and adjugate correct sign
# Now J = LDET = LA + (N - 1) * (LA - 1) = N * (LA - 1) + 1
# Hence J - LA + 1 = (N - 1) * (LA - 1) + 1 = LADJ
    ADJ = SCALE(float(2 * MOD(N, 2) - 1), N * N * (J - LA + 1), ADJ)
    for L in range(J): 
#       NL = N + L * N
        NL = N - 1 + L * N
        DET[L] = P[NL] * float(2 * MOD(N, 2) - 1)
    return (ADJ, P, DET, S)
#_______________________________________________________________________________
#_______________________________________________________________________________
def FACTOR(M, N, BS): 
    """
    FACTOR
    
    p. 197
    
    not tested yet
    """
    NX = 49
    XCOF = numpy.empty((NX + 1, ))
    COF = numpy.empty((NX + 1, ))
    RZR = numpy.empty((NX, ))    
    CZR = numpy.empty((NX, ))
    N1 = N - 1
    NZB = N1 * M
    LR = N + N1
    LAJ = (M - 1) * (LR - 1) + 1
    LDETR = (LR - 1) * M + 1
    NZR = 2 * NZB
    LRR = LAJ + N1
    R = numpy.empty((M * M * LR, ), "complex")
    AJ = numpy.empty((M * M * LAJ), "complex")
    DETR = numpy.empty((LDETR, ), "complex")
    ZR = numpy.empty((NZR, ), "complex")
    ZB = numpy.empty((M * (N - 1), ), "complex")
    RR = numpy.empty((M * M * LRR, ), "complex")
    FACT = numpy.empty((M * M * 2, ), "complex")
    DIAG = numpy.empty((M * M, ), "complex")
    CR = numpy.empty((M * M, ), "complex")
    P = numpy.empty((M * M), "complex")
    PINV = numpy.empty((M * M, ), "complex")
    V = numpy.empty((M * M, ), "complex")
    TEMP = numpy.empty((M * M * LDETR, ), "complex")
    B = numpy.empty((M * M * N, ), "complex")
    BS = numpy.empty((M * M * N, ), "complex")
    BZ1 = numpy.empty((M * M, ), "complex")
    RZ1 = numpy.empty((M * M, ), "complex")
    W = numpy.empty((M * M, ), "complex")
    BZINV = numpy.empty((M * M, ), "complex")
# complex DETP, DET, ZER, Z1
# Compute one side of autocorrelation
    HEAT(M, M, N, BS, M, M, N, BS, N, R[0, 0, N])
# Generate other side of autocorrelation
    for I in range(1, N): 
        JP = N + I - 1
        JM = N - I + 1
        for J in range(M): 
            for K in range(M): 
                JKJM = J + K * M + JM * M * M
                KJJP = K + J * M + JP * M * M
                R[JKJM] = R[KJJP]
# Find adjugate and autocorrelation of R
    POMAIN(M, LR, R, AJ, RR, DETR, TEMP)
# Find zeros of determinant
    for I in range(LDETR): 
        XCOF[I] = DETR[I].real
    POLRT(XCOF, COF, LDETR - 1, RZR, CZR, IER)
    for I in range(NZR): 
        ZR[I] = complex(RZR[I], CZR[I])
    J = 0
    for I in range(NZR): 
# It is assumed that no zero has magnitude unity
        if (abs(ZR[I]) >= 1.0): 
# CABS is abs for complex 
# Choose the zeros with magnitude greater than unity to form B
            DETR[J] = ZR[I]
            J += 1
    MOVE(NZB, DETR, ZB)
# Set B = I, FACT(J, K, 1) = 1        
# ZERO(M * M, FACT)
# ZERO(M * M, B)
    FACT = numpy.zeros((M * M, ))
    B = numpy.zeros((M * M, ))
    for J in range(M):
        JJ0 = J + J * M 
        FACT[JJ0] = 1. 
        B[JJ0] = 1. 
    LBT = 1
# Set RR = AJ
    MOVE(LAJ * M * M, AJ, RR)
    LRRT = LAJ
# Loop on the factors (I - Z * PINV * DIAG * P)
    for IFACTS in range(N1): 
# Form diagonal matrix
# ZERO(M * M, DIAG)
        DIAG = numpy.zeros((M * M, ))
        for I in range(M):
            II = I + I * M
            IIFACTS = I + IFACTS * M
            DIAG[II] = 1. / ZB[IIFACTS]
# Insert zeros in RR to get eigenvectors
        for IVECT in range(M): 
# ZERO(M * M, CR)
            CR = numpy.zeros((M * M, ))
            IVECTIFACTS = IVECT + IFACTS * M
            ZER = ZB[IVECTIFACTS]
            for I in range(LRRT):
                ZI = ZER ** (I - LAJ / 2.)
                for K in range(M):
                    for J in range(M): 
                        JK = J + K * M
                        JKI = J + K * M + I * M * M
                        CR[JK] += RR[JKI] * ZI 
            for I in range(M):
                IVECTI = IVECT + I * M
# [0, I] becomes [I]
                P[IVECTI] = CR[I]
# Form PINV * DIAG * P
# [1, 1, 2] becomes [1]
        BRAINY(M, M, 1, DIAG, M, M, 1, P, FACT[1])
        FADDEJ(M, P, PINV, DETP, DIAG, TEMP)
# [1, 1, 2] becomes [1]
        MOVE(M * M, FACT[1])
# [1, 1, 2] becomes [1]
        BRAINY(M, M, 1, PINV, M, M, 1, P, FACT[1])
# [1, 1, 2] becomes [1]
        SCALE(-1., M * M, FACT[1])
# Set RR = RR * (I - Z * PINV * DIAG * P)
        BRAINY(M, M, LRRT, RR, M, M, 2, FACT, TEMP)
        LRRT += 1
        MOVE(LRRT * M * M, TEMP, RR)
# Set B = B * (I - FACT)
        BRAINY(M, M, LBT, B, M, M, 2, FACT, TEMP)
        LBT += 1
        MOVE(M * M * LBT, TEMP, B)
    HEAT(M, M, N, B, M, M, N, B, N, TEMP)
# Form B(Z = 1) and R(Z = 1)
    RZ1 = numpy.zeros((M * M, ))
    BZ1 = numpy.zeros((M * M, ))
    for J in range(M): 
        for K in range(M): 
            JK = J + K * M
            for I in range(N):
                JKI = JK + I * M
                BZ1[JK] += B[JKI]
            for I in range(LR): 
                JKI = JK + I * M
                RZ1[JK] += R[JKI]
# Form BINV(Z = 1) * R(Z = 1) * BINVTRANSP(Z = 1) = W
    FADDEJ(M, BZ1, BZINV, DET, TEMP, W)
    BRAINY(M, M, 1, BZINV, M, M, 1, RZ1, TEMP)
    HEAT(M, M, 1, TEMP, M, M, 1, BZINV, 1, W)
# Triangularize W
    TRIANG(M, W, V, TEMP)
# Put B in causal-chain form
    HEAT(M, M, N, B, M, M, 1, V, N, TEMP)
    MOVE(M * M * N, TEMP, B)
    HEAT(M, M, N, B, M, M, N, B, N, TEMP)
#    return (LR, R, LAJ, AJ, LDETR, DETR, NZR, ZR, ZB, LRR, RR, FACT, DIAGCR, P, PINV, V, TEMP, B, BZ1, RZ1, W, BZINV)
    return (LR, R, LAJ, AJ, LDETR, DETR, NZR, ZR, ZB, B)
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
    return (C, LC)
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


