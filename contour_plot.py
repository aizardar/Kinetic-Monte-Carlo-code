import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import Rbf
from matplotlib import cm
import scipy.ndimage as ndimage


#Z = np.cos((X**2+Y**2)/200.)+ np.random.normal(size=X.shape)

# Increase the value of sigma to increase the amount of blurring.
# order=0 means gaussian kernel




x, y, z = np.genfromtxt(r'D_vs_Al_concentration_250K.dat', unpack=True)
logz = np.log10(z)

tix = np.linspace(x.min(), x.max(), 48)
tiy = np.linspace(y.min(), y.max(), 48)
xi,yi = np.meshgrid(tix,tiy)
Z = interpolate.griddata((x, y), logz, (xi[None,:], yi[:,None]), method='cubic')
Z2 = ndimage.gaussian_filter(Z, sigma=1.0, order=0)
#rbf = Rbf(x, y, logz, epsilon=10)
#ZI = rbf(XI, YI)
plt.figure()
fig=plt.figure()
ax=fig.add_subplot(1,2,1)
ax.imshow()
ax=fig.add_subplot(1,2,2)
ax.imshow(Z2)

#plt.pcolor(XI, YI, ZI, cmap=cm.jet)


#cp = plt.contourf(xi, yi, zi)
#plt.pcolor(xi,yi,zi)
#plt.colorbar()
#plt.xlabel("Al Concentration")
#plt.ylabel("Mn Concentration")
#plt.savefig('1000')
plt.show()

