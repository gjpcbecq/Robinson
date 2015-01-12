import numpy
SIN = numpy.sin
COS = numpy.cos

#_______________________________________________________________________________
def TRIG(LX, X, W): 
    """
TRIG computes one value of Fourier transform by the sum of angles TRIGonometric formula for sine and cosine. 
"""
    COSNW = 1.
    SINNW = 1.
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
    return (C, LC)
#_______________________________________________________________________________
#_______________________________________________________________________________
def POMAIN(N, LA, A, ADJ, P, DET, S): 
    """
    POMAIN Polynomial matrix inversion
    
    p. 162
    """
#    complex A, ADJ, P, DET, S
    MOVE(N * N * LA, A, S)
    J = LA
    for L in range(N): 
# Calculate coefficients P[., K] of characteristic polynomial
        for K in range(J):
            LK = L + K * N
            P[LK] = 0.
            for I in range(N): 
                IIK = I + I * N + K * N * N
                P[LK] += S[IIK] / float(L)
        if (L != N): 
# Substract P[., K]*identity matrix 
            MOVE(N * N * J, S, ADJ)
            for I in range(N): 
                for K in range(J): 
                    IIK = I + I * N + K * N * N
                    LK = L + K * N
                    ADJ[IIK] -= P[LK]
# Multiply by input matrix 
            BRAINY(N, N, LA, A, N, N, J, ADJ, S)
            #J += LA - 1
            J += LA
# Give determinant and adjugate correct sign
        SCALE(float(2 * MOD(N, 2)), N * N * (J - LA + 1), ADJ)
        for L in range(J): 
            NL = N + L * N
            DET[L] = P[NL] * float(2 * MOD(N, 2) - 1)
    return (ADJ, P, DET, S)
    
