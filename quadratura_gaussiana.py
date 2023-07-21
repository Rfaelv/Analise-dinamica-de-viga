import numpy as np


def quadraturaGaussiana(n, inf, sup, _x, _y):
    pontos, pesos = np.polynomial.legendre.leggauss(n)
    sol = 0

    for i in range(n):
        x  = (sup - inf) / 2 * pontos[i] + (sup + inf) / 2
        yi = _y[np.absolute(_x-x).argmin()-1]
        yf = _y[np.absolute(_x-x).argmin()+1]
        xi = _x[np.absolute(_x-x).argmin()-1]
        xf = _x[np.absolute(_x-x).argmin()+1]
        y  = _y[np.absolute(_x-x).argmin()-1]+ (yf-yi)/(xf-xi)*(x-xi) 

        value = (sup - inf) / 2 * y
        sol  += pesos[i] * value

    return sol