from scipy.fftpack import fft2
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math



gridParameters = np.genfromtxt('gridParameters.txt')
lengthOfSimulationBoxInX = gridParameters[6]
lengthOfSimulationBoxInY = gridParameters[7]
lengthOfOneBoxInX = gridParameters[0] * gridParameters[3]
lengthOfOneBoxInY = gridParameters[1] * gridParameters[4]
numberOfBoxesInX =  lengthOfSimulationBoxInX / lengthOfOneBoxInX
numberOfBoxesInY =  lengthOfSimulationBoxInY / lengthOfOneBoxInY
numberOfGridPointsInX = numberOfBoxesInX * gridParameters[3]

t = np.linspace(0,lengthOfSimulationBoxInX,numberOfGridPointsInX)
dt = t[1] - t[0]
sampleFrequency = 1.0/dt
nyquistFrequency = sampleFrequency / 2.0

Eges = np.genfromtxt('E_fields/E_field0.txt') #signal
Sx = np.genfromtxt('FourierAnalysis/Sx.txt') #signal
Sy = np.genfromtxt('FourierAnalysis/Sy.txt') #signal
Sz = np.genfromtxt('FourierAnalysis/Sz.txt') #signal
Yx = np.fft.fft2(Sx)
Yy = np.fft.fft2(Sy)
Yz = np.fft.fft2(Sz)
N = len(Yx)/2+1
# SxMax = 2 * pow(10,-13)
# SyMax = 1 * pow(10,-12)
# SzMax = 5 * pow(10,-15)
SxMax = Sx.max()
SyMax = Sy.max()
SzMax = Sz.max()

scaleFactor = pow(10,-8)
fig = plt.figure(figsize=(12, 6))
# plt.subplot(231)
# axp = plt.imshow(Eges, aspect='equal', origin='lower', cmap = 'jet', extent=(0,lengthOfSimulationBoxInX,0,lengthOfSimulationBoxInY), vmin=0, vmax=colorbarLimit)
# plt.colorbar()
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title('$|E|^2$')

plt.subplot(231)
axp = plt.imshow(Sx, origin='lower', cmap = 'jet', extent=(0,lengthOfSimulationBoxInX,0,lengthOfSimulationBoxInY), vmin=-SxMax, vmax=SxMax)
plt.colorbar(format='%.0e')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('$S_x$')

plt.subplot(232)
axp = plt.imshow(Sy, origin='lower', cmap = 'jet', extent=(0,lengthOfSimulationBoxInX,0,lengthOfSimulationBoxInY), vmin=-SyMax, vmax=SyMax)
plt.colorbar(format='%.0e')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('$S_y$')

plt.subplot(233)
axp = plt.imshow(Sz, origin='lower', cmap = 'jet', extent=(0,lengthOfSimulationBoxInX,0,lengthOfSimulationBoxInY), vmin=-SzMax, vmax=SzMax)
plt.colorbar(format='%.0e')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('$S_z$')

plt.subplot(234)
plt.imshow(2 * np.abs(Yx[:N,:N]) / N, origin='lower', cmap = 'jet', extent=(0,nyquistFrequency,0,nyquistFrequency))
plt.colorbar(format='%.0e')
plt.xlabel('$\omega_x$')
plt.ylabel('$\omega_y$')
plt.title('FFT $S_x$')
plt.xlim([0,2 * scaleFactor])
plt.ylim([0,2 * scaleFactor])

plt.subplot(235)
plt.imshow(2 * np.abs(Yy[:N,:N]) / N, origin='lower', cmap = 'jet', extent=(0,nyquistFrequency,0,nyquistFrequency))
plt.colorbar(format='%.0e')
plt.xlabel('$\omega_x$')
plt.ylabel('$\omega_y$')
plt.title('FFT $S_y$')
plt.xlim([0,2 * scaleFactor])
plt.ylim([0,2 * scaleFactor])

plt.subplot(236)
plt.imshow(2 * np.abs(Yz[:N,:N]) / N, origin='lower', cmap = 'jet', extent=(0,nyquistFrequency,0,nyquistFrequency))
plt.colorbar(format='%.0e')
plt.xlabel('$\omega_x$')
plt.ylabel('$\omega_y$')
plt.title('FFT $S_z$')
plt.xlim([0,2 * scaleFactor])
plt.ylim([0,2 * scaleFactor])

#plt.subplots_adjust(wspace=0.3, hspace=0.3)
fig.tight_layout()

#plt.show()

# ax1.imshow(y, aspect='auto', origin='lower', cmap = 'jet', extent=(0,32,0,32))
# ax2.imshow(np.real(Z), aspect='auto', origin='lower', cmap = 'jet', extent=(0,32,0,32))
filename = 'fourierAnalysisPoynting2D'
fig.savefig("{}.png".format(filename),bbox_inches='tight',dpi=300)
