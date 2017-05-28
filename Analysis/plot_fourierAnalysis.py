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
Ex = np.genfromtxt('FourierAnalysis/Ex/E_field_x0.txt') #signal
Ey = np.genfromtxt('FourierAnalysis/Ey/E_field_y0.txt') #signal
Ez = np.genfromtxt('FourierAnalysis/Ez/E_field_z0.txt') #signal
Yx = np.fft.fft2(Ex)
Yy = np.fft.fft2(Ey)
Yz = np.fft.fft2(Ez)
N = len(Yx)/2+1
# ExMax = 8 * pow(10,-7)
# EyMax = 1 * pow(10,-6)
# EzMax = 4 * pow(10,-9)
ExMax = Ex.max()
EyMax = Ey.max()
EzMax = Ez.max()

scaleFactor = pow(10,-8)
fig = plt.figure(figsize=(12, 6))
# plt.subplot(231)
# axp = plt.imshow(Eges, aspect='equal', origin='lower', cmap = 'jet', extent=(0,lengthOfSimulationBoxInX,0,lengthOfSimulationBoxInY), vmin=0, vmax=colorbarLimit)
# plt.colorbar()
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title('$|E|^2$')

plt.subplot(231)
axp = plt.imshow(Ex, origin='lower', cmap = 'jet', extent=(0,lengthOfSimulationBoxInX,0,lengthOfSimulationBoxInY), vmin=-ExMax, vmax=ExMax)
plt.colorbar(format='%.0e')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('$E_x$')

plt.subplot(232)
axp = plt.imshow(Ey, origin='lower', cmap = 'jet', extent=(0,lengthOfSimulationBoxInX,0,lengthOfSimulationBoxInY), vmin=-EyMax, vmax=EyMax)
plt.colorbar(format='%.0e')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('$E_y$')

plt.subplot(233)
axp = plt.imshow(Ez, origin='lower', cmap = 'jet', extent=(0,lengthOfSimulationBoxInX,0,lengthOfSimulationBoxInY), vmin=-EzMax, vmax=EzMax)
plt.colorbar(format='%.0e')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('$E_z$')

plt.subplot(234)
plt.imshow(2 * np.abs(Yx[:N,:N]) / N, origin='lower', cmap = 'jet', extent=(0,nyquistFrequency,0,nyquistFrequency))
plt.colorbar(format='%.0e')
plt.xlabel('$\omega_x$')
plt.ylabel('$\omega_y$')
plt.title('FFT $E_x$')
plt.xlim([0,2 * scaleFactor])
plt.ylim([0,2 * scaleFactor])

plt.subplot(235)
plt.imshow(2 * np.abs(Yy[:N,:N]) / N, origin='lower', cmap = 'jet', extent=(0,nyquistFrequency,0,nyquistFrequency))
plt.colorbar(format='%.0e')
plt.xlabel('$\omega_x$')
plt.ylabel('$\omega_y$')
plt.title('FFT $E_y$')
plt.xlim([0,2 * scaleFactor])
plt.ylim([0,2 * scaleFactor])

plt.subplot(236)
plt.imshow(2 * np.abs(Yz[:N,:N]) / N, origin='lower', cmap = 'jet', extent=(0,nyquistFrequency,0,nyquistFrequency))
plt.colorbar(format='%.0e')
plt.xlabel('$\omega_x$')
plt.ylabel('$\omega_y$')
plt.title('FFT $E_z$')
plt.xlim([0,2 * scaleFactor])
plt.ylim([0,2 * scaleFactor])

#plt.subplots_adjust(wspace=0.3, hspace=0.3)
fig.tight_layout()

#plt.show()

# ax1.imshow(y, aspect='auto', origin='lower', cmap = 'jet', extent=(0,32,0,32))
# ax2.imshow(np.real(Z), aspect='auto', origin='lower', cmap = 'jet', extent=(0,32,0,32))
filename = 'fourierAnalysis2D'
fig.savefig("{}.png".format(filename),bbox_inches='tight',dpi=300)
