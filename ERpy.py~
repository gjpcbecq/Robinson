"""
ERpy binding to ER
"""
#_______________________________________________________________________________
import ER
import numpy
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

