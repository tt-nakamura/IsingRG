import numpy as np
import matplotlib.pyplot as plt
from IsingRG import zoom_K

N = 50
K1,K2 = 0.001, 10

K = zoom_K(K1,N)
plt.plot(K[:-1], K[1:], 'r.', label=r'$K_0=%g$'%K1)
print(K[-1])

K = zoom_K(K2,N)
plt.plot(K[:-1], K[1:], 'b.', label=r'$K_0=%g$'%K2)
print(K[-1])

plt.plot([K1,K2],[K1,K2], ':', label=r'$K_{n+1}=K_n$')
plt.axis([0.504, 0.51, 0.504, 0.51])
plt.xlabel(r'$K_n$', fontsize=14)
plt.ylabel(r'$K_{n+1}$', fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
