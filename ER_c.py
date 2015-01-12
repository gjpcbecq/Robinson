"""
Enders A. Robinson with C functions 

"""
import numpy
import ctypes
import os
path = __file__.rpartition(os.sep)[0]
print(path)
ffnlib = os.path.join(path, "ER.dylib")
print(ffnlib)
libc = ctypes.CDLL(ffnlib)
print(libc)
print(dir(libc))
NDPOINTERDOUBLE = numpy.ctypeslib.ndpointer(dtype=numpy.double)
libc.QUADCO.argtypes = [ctypes.c_int, ctypes.c_int, NDPOINTERDOUBLE, NDPOINTERDOUBLE, NDPOINTERDOUBLE]
libc.COHERE.argtypes = [ctypes.c_int, ctypes.c_int, NDPOINTERDOUBLE, NDPOINTERDOUBLE]
libc.MACRO.argtypes = [ctypes.c_int, ctypes.c_int, NDPOINTERDOUBLE, ctypes.c_int, NDPOINTERDOUBLE, ctypes.c_int, NDPOINTERDOUBLE]
libc.REMAV.argtypes = [ctypes.c_int, NDPOINTERDOUBLE]

def MACRO(N, LX, X, LY, Y, LG, G): 
    """
    G = numpy.empty((LG*N*N, ), dtype="float")
    """
    libc.MACRO(N, LX, X, LY, Y, LG, G)
    return G

def QUADCO(L, N, R, S, SP):
    """
    S = numpy.empty((L*N*N, ), dtype="float")
    SP = numpy.empty((L*N*N, ), dtype="float")
    """
    libc.QUADCO(L, N, R, S, SP)
    return S
    
def COHERE(L, N, S, C):
    """C = numpy.empty((L*N*N), dtype="float")
    """
    libc.COHERE(L, N, S, C)
    return C

def REMAV(LY, Y):
    libc.REMAV(LY, Y)
    return Y


