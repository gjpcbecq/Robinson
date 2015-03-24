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

