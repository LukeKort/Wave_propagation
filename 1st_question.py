# Cálculo da transformada de Fourier

import numpy as np
import matplotlib.pyplot as plt

rho = 7850  # Bar's density
E = 210e9 # Bar's Young Modulos
c_0 = np.sqrt(E/rho) # Bar's propagation velocity

N = 10 # Number of discratization points

A = 1 # Maxium amplitude

a = -2 # Sign limits
b = 0
c = 1

x_disc = np.linspace(-2,1,N) # Sign limits

y_disc = np.piecewise(x_disc,
                      [x_disc < a, 
                       (x_disc >= a) & (x_disc <= b), 
                       (x_disc > b) & (x_disc <= c),
                       (x_disc > c)
                       ],
                      [0,
                       lambda x: A*(x/2 + 1),
                       lambda x: A*(-x + 1),
                       0
                       ]
                      ) # Discretizes the sign along the x axis

# Process and recovery sign

from scipy.fft import fft, fftfreq, ifft

# Processing of the FFT

yf = fft(y_disc) # Calculates the amplitudes
xf = fftfreq(N) # Calculates the frequencies associated to the amplitudes above

plt.figure('Frequency domain')
plt.scatter(xf, np.abs(yf))
plt.xlabel('Número de onda')
plt.ylabel('Amplitude')
plt.legend()

# Recovery of the sign

recoverd_sign = ifft(yf) # Recuperates the original sign from the FFT frequencies

plt.figure('Original and Recovered Sign')
plt.plot(x_disc, y_disc, label = 'Original Sign', linestyle = '-') # Original Sign
plt.plot(x_disc, recoverd_sign, label = 'Recovered sign', linestyle = 'dotted') # Recovered sign
plt.xlabel('x')
plt.ylabel('u(x,t)')
plt.legend()

# Calculate the wave propagation

A = 0.5*yf

k = 2*np.pi*xf

omega = c_0*k

x_linspace = np.linspace(-10, 10, 10000)

u_x_t = np.zeros(np.shape(x_linspace)[0])

def u(x,t):
    u = 0

    i = int(N/2)

    for n in range(i):
        u =+ (yf[n])*np.exp(1j*(k[n]*x - omega[n]*t)) + (yf[n + i])*np.exp(-1j*(k[n + i]*x + omega[n + i]*t))

    return np.abs(u)

for i in range(np.shape(x_linspace)[0]):
    u_x_t[i] = u(x_linspace[i],1)

plt.figure('Wave propagation')
plt.plot(x_linspace, u_x_t)
plt.show()