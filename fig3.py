import numpy as np
import matplotlib.pyplot as plt
from IsingRG import zoom_KUC

N = 100
K01,U01,C01 = 0.001, 0, 0
K02,U02,C02 = 10, -2, 0

plt.figure(figsize=(6.4,6))

K1,U1,C1 = zoom_KUC(K01,U01,C01,N)
K2,U2,C2 = zoom_KUC(K02,U02,C02,N)

plt.subplot(211)
plt.axis([0.504, 0.51, -2.33, -2.2])
plt.ylabel(r'$u$ = internal enegy', fontsize=14)
plt.plot(K1, U1, 'r.', label=r'$K_0=%g$'%K01)
plt.plot(K2, U2, 'b.', label=r'$K_0=%g$'%K02)
plt.legend(fontsize=14)

plt.subplot(212)
plt.axis([0.504, 0.51, 0, 50])
plt.ylabel(r'$c$ = specific heat', fontsize=14)
plt.plot(K1, C1, 'r.', label=r'$K_0=%g$'%K01)
plt.plot(K2, C2, 'b.', label=r'$K_0=%g$'%K02)
plt.legend(fontsize=14)
plt.xlabel(r'$K = J/(k_{\rm B}T)$', fontsize=14)

plt.tight_layout()
plt.savefig('fig3.eps')
plt.show()
