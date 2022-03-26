import numpy as np
import matplotlib.pyplot as plt

fname = ['fig4_K3.txt', 'fig4_K3z.txt',
         'fig4_Kc.txt', 'fig4_Kcz.txt',
         'fig4_K5.txt', 'fig4_K5z.txt']
title = [r'$K = 0.3$', 'zoomed out',
         r'$K = 0.4407$', 'zoomed out',
         r'$K = 0.5$', 'zoomed out']

plt.figure(figsize=(4,6))

for i,fname in enumerate(fname):
    s = np.loadtxt(fname)
    plt.subplot(3,2,i+1)
    plt.title(title[i])
    plt.imshow(s, cmap='gray')
    plt.xticks([])
    plt.yticks([])

plt.tight_layout()
plt.savefig('fig4.eps')
plt.show()
