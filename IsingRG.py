# referece: H. J. Maris and L. P. Kadanoff
#   American Journal of Physics 46 (1978) 652

import numpy as np
from scipy.optimize import newton
from scipy.special import ellipk,ellipe

# K = J/(kT) (T = temperature, J = spin interaction)
# critical value of K
Kc = newton(lambda x: x - 3/8*np.log(np.cosh(4*x)), 1,
                lambda x: 1 - 3/2*np.tanh(4*x))

Kc_exact = np.log(1 + np.sqrt(2))/2 # Onsager

# critical exponent of specific heat (|K-Kc|^-alpha)
alpha = 2 - 1/np.log2(1.5*np.tanh(4*Kc))

def zoom_K(K_init, N):
    """ renormalization group iteration of two dimensional
    Ising model, converging to fixed point
    K_init = initial value of K = J/(kT)
    N = number of iterations
    return K[0:N]
    """
    K1,K = K_init,[K_init]
    for _ in range(N):
        K1 = np.arccosh(np.exp(8*K1/3))/4
        K.append(K1)
    return np.asarray(K)

def zoom_KUC(K_init, U_init, C_init, N):
    """ renormalization group iteration of two dimensional
    Ising model, converging to fixed point
    K_init, U_init, C_init = initial values of
      K = J/(kT),
      U = internal energy per spin / J,
      C = specific heat per spin / k
    N = number of iterations
    return K[0:N],U[0:N],C[0:N]
    """
    def KUC_next(K,U,C):# one iteration
        cK = np.exp(8*K/3)
        K1 = np.arccosh(cK)/4
        tK = np.tanh(4*K1)
        U1 = (3*U - 3 + 2/(cK+1))*tK/4
        C1 = ((9/8*C/K**2 + 2*cK/(cK+1)**2)*tK**2
              - 4*U1/tK/cK**2)*K1**2
        return K1,U1,C1

    K,U,C = K_init, U_init, C_init
    y = [[K,U,C]]
    for _ in range(N):
        K,U,C = KUC_next(K,U,C)
        y.append([K,U,C])
    return np.asarray(y).T

def exact_UC(K):
    """ Onsager's exact solution of
    two dimensional Ising model
    K = J/(kT)
    return U,C (same shape as K)
      U = internal energy per spin / J
      C = specific heat per spin / k
    """
    tK = np.tanh(2*K)
    k = 2*tK/np.cosh(2*K)
    k1 = 2*tK**2 - 1
    Kk = ellipk(k)*2/np.pi
    Ek = ellipe(k)*2/np.pi
    U = -(1 + k1*Kk)/tK
    C = (K/tK)**2*(2*(Kk - Ek) - (1-k1)*(4/np.pi**2 + k1*Kk))
    return U,C
