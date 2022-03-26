import numpy as np
import matplotlib.pyplot as plt
from IsingRG import exact_UC

K,U1,C1 = np.loadtxt('fig5_L32.txt').T
K,U2,C2 = np.loadtxt('fig5_L128.txt').T
Ke = np.linspace(K[0], K[-1], 256)
Ue,Ce = exact_UC(Ke)

plt.figure(figsize=(6.4,6))
plt.subplot(211)
plt.plot(K,U1,'.-', label=r'$L=32$')
plt.plot(K,U2,'.-', label=r'$L=128$')
plt.plot(Ke,Ue, label='exact')
plt.ylabel(r'$u$ = internal energy', fontsize=14)
plt.legend(fontsize=14)

plt.subplot(212)
plt.plot(K,C1,'.-', label=r'$L=32$')
plt.plot(K,C2,'.-', label=r'$L=128$')
plt.plot(Ke,Ce, label='exact')
plt.ylabel(r'$c$ = specific heat', fontsize=14)
plt.xlabel(r'$K = J/(k_{\rm B}T)$', fontsize=14)
plt.legend(fontsize=14)

plt.tight_layout()
plt.savefig('fig5.eps')
plt.show()
