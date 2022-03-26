import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(2,2))
x,y = np.mgrid[1:3:2j, 0:2:2j]; plt.plot(x,y,'ko')
x,y = np.mgrid[0:2:2j, 1:3:2j]; plt.plot(x,y,'ko')
x,y = np.mgrid[0:2:2j, 0:2:2j]; plt.plot(x,y,'ko', fillstyle='none')
x,y = np.mgrid[1:3:2j, 1:3:2j]; plt.plot(x,y,'ko', fillstyle='none')
plt.text(1, 2.9, '1', va='top', ha='center')
plt.text(0, 1.9, '3', va='top', ha='center')
plt.text(1, 0.9, '5', va='top', ha='center')
plt.text(2, 1.9, '7', va='top', ha='center')
plt.text(1, 1.9, '2', va='top', ha='center')
plt.text(0, 0.9, '4', va='top', ha='center')
plt.text(1,-0.1, '6', va='top', ha='center')
plt.text(2, 0.9, '8', va='top', ha='center')
plt.box('off')
plt.xticks([])
plt.yticks([])
plt.axis('equal')
plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
