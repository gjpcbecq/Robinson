"""
ERpy binding to ER
"""
#_______________________________________________________________________________
import ER
import numpy
zeros = numpy.zeros
import pylab
import ER_c
#_______________________________________________________________________________
#_______________________________________________________________________________
def POLRT(A): 
	"""
	POLRT
	
	p.35
	"""
	XCOF = A
	COF = XCOF.copy()
	M = XCOF.size - 1
	ROOTR = numpy.zeros(M)
	ROOTI = numpy.zeros(M)
	IER = 0
	
	(COF, ROOTR, ROOTI, IER) = ER.POLRT(XCOF, COF, M, ROOTR, ROOTI, IER)
	return (ROOTR, ROOTI, IER)
#_______________________________________________________________________________
#_______________________________________________________________________________
def POMAIN(N, LA, A): 
	"""
	POMAIN
	
	p. 162
	"""
	(ADJ, P, DET, S) = ER.POMAIN(N, LA, A, [], [], [], [])
	
	return (ADJ, P, DET)
#_______________________________________________________________________________
#_______________________________________________________________________________
def POMAEVAL(N, LA, A, z0): 
    """
    eval polynomial matrix at z0
    """
    z0 = numpy.complex(z0)
    P = numpy.zeros(N * N, "complex")
    for I in range(N): 
        for J in range(N): 
            for K in range(LA):
                ZK = z0 ** K
                IJ = I + J * N
                IJK = I + J * N + K * N * N
                P[IJ] += A[IJK] * ZK
    return P
#_______________________________________________________________________________
#_______________________________________________________________________________
def WIENER(N, LX, X, M, LZ, Z, LR, LW, FLOOR): 
    """
    WIENER filter 
    
    (F, E, Y) = WIENER(N, LX, X, M, LZ, Z, LR, LW, FLOOR)
    
     N: number of input channels to filter >= 1
    LX: length of input time series
     X: N-channel desired input time series
         array 1d in multiplexed mode:  
         [x1(1), x2(1), ..., xN(1), x1(2), x2(2), ..., xN(2), ..., xN(LX)] 
     M: number of output channels from filter Z >= 1
    LZ: length of desired output time series 
     Z: M-channel desired output time series
    
    """
    LF = M * N * LR
    F = numpy.array(LF)
    E = 0
    LY = 10
    Y = numpy.array(LY)
    LS = N * N * (5 * LR + 6) + M * N * (LR + 2) + 2 * M * M
    S = numpy.array(LS)
    (LF, F, E, LY, Y, S) = ER.WIENER_1(N, LX, X, M, LZ, Z, LR, LW, FLOOR, LF, F, E, LY, Y, S)
    
    return (F, Y, E)
#_______________________________________________________________________________
#_______________________________________________________________________________
def getTraceMode(N, LX, X): 
    """
    get trace mode from X in multiplexed mode 
    
    Y = getTraceMode(N, LX, X)

"""
    Y = numpy.empty(N * LX)
    for I in range(N): 
        for J in range(LX):
            IJ = I + J * N
            JI = J + I * LX
            Y[JI] = X[IJ]
    return Y

MMTOTM = getTraceMode
#_______________________________________________________________________________
#_______________________________________________________________________________
def getMultiplexedMode(N, LX, X): 
    """
    get multiplexed mode from X in trace mode 

    Y = getMultiplexedMode(N, LX, X)

    """
    Y = numpy.empty(N * LX)
    for I in range(N): 
        for J in range(LX):
            IJ = I + J * N
            JI = J + I * LX
            Y[IJ] = X[JI]
    return Y

TMTOMM = getMultiplexedMode
#_______________________________________________________________________________
#_______________________________________________________________________________
def NDTOMM(X): 
    """
    nd array to multiplexed mode
    
    X shape is (N, LX)
    
    """
    (N, LX) = X.shape
    Y = X.flatten()
    Z = getMultiplexedMode(N, LX, Y)
    return (N, LX, Z)    
#_______________________________________________________________________________
#_______________________________________________________________________________
def NDTOTM(X): 
    """
    nd array to trace mode
    
    X shape is (N, LX)
    """
    (N, LX) = X.shape
    Y = X.flatten()
    return (N, LX, Y)    
#_______________________________________________________________________________
#_______________________________________________________________________________
def TMTOND(N, LX, X): 
    """
    trace mode to nd array
    
    """
    Y = X.reshape(N, LX)
    return Y
#_______________________________________________________________________________
#_______________________________________________________________________________
def MMTOND(N, LX, X): 
    """
    multiplexed mode to nd array
    
    """
    Y = getTraceMode(N, LX, X)
    Z = Y.reshape(N, LX)
    return Z
#_______________________________________________________________________________
#_______________________________________________________________________________
def SPIKER(B, LA): 
    """
    """
    LB = B.size
    A = numpy.zeros(LA)
    LC = LB + LA - 1
    C = numpy.zeros(LC)
    INDEX = 0
    ERRORS = numpy.zeros(LC)
    SPACE = numpy.zeros(3 * LA)
    (A, LC, C, INDEX, ERRORS) = ER.SPIKER(LB, B, LA, A, LC, C, INDEX, ERRORS, SPACE)
    return (A, C, INDEX, ERRORS)
#_______________________________________________________________________________
#_______________________________________________________________________________
def SHAPER(B, D, LA): 
    """
    """
    LB = B.size
    LD = D.size
    A = numpy.zeros(LA)
    LC = LB + LA - 1
    LCD = LC + LD - 1
    C = numpy.zeros(LCD)
    INDEX = 0
    ERRORS = numpy.zeros(LCD)
    SPACE = numpy.zeros(3 * LA)
    (A, LC, C, INDEX, ERRORS, S) = ER.SHAPER(LB, B, LD, D, LA, A, LC, C, INDEX, ERRORS, SPACE)
    return (A, C, INDEX, ERRORS)
#_______________________________________________________________________________
#_______________________________________________________________________________
def SHAPE(B, D, LA): 
    """
    """
    LB = B.size
    LD = D.size
    A = numpy.zeros(LA)
    LC = LB + LA - 1
    LCD = LC + LD - 1
    C = numpy.zeros(LCD)
    INDEX = 0
    ASE = numpy.zeros(LCD)
    SPACE = numpy.zeros(3 * LA)
    (A, LC, C, ASE, SPACE) = ER.SHAPE(LB, B, LD, D, LA, A, LC, C, ASE,  SPACE)
    return (A, C, ASE)
#_______________________________________________________________________________
#_______________________________________________________________________________
def MACRO(X, Y, LG): 
    """
    MACRO multichannel cross correlation
    
    (G, N) = MACRO(X, Y, LG) 
        
    X: (nDimX, nObsX) (N, LX)
    Y: (nDimY, nObsY) (N, LY)
    
    Output
    (G, N) 
    
    """
    (N, LX, X) = NDTOTM(X)
    (N, LY, Y) = NDTOTM(Y)
    # print(X)
    # pylab.plot(X)
    G = zeros((LG * N * N))
    # G = ER.MACRO(N, LX, X, LY, Y, LG, G)
    G = ER_c.MACRO(N, LX, X, LY, Y, LG, G)
    return (G, N)
#_______________________________________________________________________________
#_______________________________________________________________________________
def MACRO_partial(N, LG, G, I, J, K): 
    """
    MACRO_partial multichannel cross partial correlation
    
    G = MACRO(X, Y, LG)
    
    Example
    
    LG = 4
    N = 3
    G = array([
        1., 0.9, 0.8, 0.7, 0.2, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 
        0.1, 0.1, 0.2, 0.1, 1, 0.8, 0.6, 0.4, 0.1, 0., 0., 0., 
        0., 0., 0., 0., 0.2, 0.3, .4, .5, 1., 0.7, 0.4, 0.1])
    GP = MACRO_partial(N, LG, G, 0, 1, 2)
    print(GP)
    """
    GP = numpy.zeros(LG) 
    for IG in range(LG): 
        IGIJ = IG + I * N + J * N * N
        IGIK = IG + I * N + K * N * N
        IGJK = IG + J * N + K * N * N
        den1 = (1. - G[IGIK]) ** 0.5
        den2 = (1. - G[IGJK]) ** 0.5 
        den = den1 * den2
        num = G[IGIJ] - G[IGIK] * G[IGJK]
        GP[IG] = num / den
    return GP
#_______________________________________________________________________________
#_______________________________________________________________________________
def Sxx(X, L): 
    """
    inter spectra
    
    $$ \Phi_{kj} (f) = C_{kj}(f) - i \, Q_{kj}(f) $$
    
    $$ \Phi_{kj} (f) = \Phi_{jk}^{*}(f) = C_{jk}(f) + i \, Q_{kj}(f) $$
    
    For power spectrum multiply the array by conjugate.
    
    """
    (N, LX, X) = NDTOTM(X)
    Rxx = zeros((L * N * N))
    Rxx = ER.MACRO(N, LX, X, LX, X, L, Rxx)
    S0 = ER.QUADCO(L, N, Rxx)
    S1 = zeros((L, N, N), 'complex')
    for i in range(N): 
        for j in range(N): 
            for k in range(L): 
                ijk = L * N * i + j * L + k
                S1[k, i, j] = S0[ijk]
    S = zeros((L, N, N), 'complex')
    for iF in range(L): 
        for i in range(N): 
            for j in range(i, N): 
                if (i == j): 
                    S[iF, i, j] = S1[iF, i, j]
                else : 
                    S[iF, i, j] = S1[iF, i, j] - complex(0, j) * S1[iF, j, i] 
                    S[iF, j, i] = S[iF, i, j].conj() 
    return (S, S1, N)    
#_______________________________________________________________________________
#_______________________________________________________________________________
def plot_MACRO(G, LG, NX, xt=[], xl=[], yt=[], rangeXY='', mode='oneSide'): 
    """
    plot_MACRO(G, LG, NX)
    
    
    """
    for i in range(NX): 
        for j in range(NX):
            JIPO = j + NX * i + 1
            pylab.subplot(NX, NX, JIPO)
            ZJI = j * LG + LG * NX * i
            ZIJ = i * LG + LG * NX * j
            Gl = G[ZIJ: ZIJ + LG][::-1]
            Gr = G[ZJI: ZJI + LG]
            if mode == 'oneSide': 
                H = pylab.r_[Gr]
            else : 
                H = pylab.r_[Gl, Gr[1:]]
            pylab.plot(H, '.-')
            pylab.xticks(xt, xl)
            pylab.yticks(yt, yt)
            if (rangeXY != ''): 
                pylab.axis(rangeXY)
#_______________________________________________________________________________
