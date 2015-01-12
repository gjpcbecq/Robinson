def eval_pol_and_deriv(COF, X, NX, NXX, Y, SUMSQ, ALPHA, N):
    # evaluate polynomial and derivatives    
    UX = 0.0
    UY = 0.0
    V = 0.0
    YT = 0.0
    XT = 1.0
    U = COF[N]
    print(U)
    if (U == 0): 
        X = 0.0
        NX -= 1
        NXX -= 1
        Y = 0.0
        SUMSQ = 0.0
        ALPHA = X
        N -= 1
        return (X, NX, NXX, Y, SUMSQ, ALPHA, N, 0, 0)
    for I in range(N):
        print((XT, YT, XT2, YT2, X, Y, U, V, UX, UY))
        L = N - I - 1
        XT2 = X * XT - Y * YT
        YT2 = X * YT + Y * XT
        U += COF[L] * XT2
        V += COF[L] * YT2
        FI = I
        UX += FI * XT * COF[L]
        UY -= FI * YT * COF[L]
        XT = XT2
        YT = YT2
    SUMSQ = UX * UX + UY * UY
    print(SUMSQ)
    if (SUMSQ != 0): 
        DX = (V * UY - U * UX) / SUMSQ
        X += DX
        DY = -(U * UY + V * UX) / SUMSQ
        Y += DY
    return (X, NX, NXX, Y, SUMSQ, ALPHA, N, DX, DY)
    
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
    