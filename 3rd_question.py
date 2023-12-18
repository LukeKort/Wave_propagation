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

# Calculate the wave propagation

A = 0
B = 0

for i in range(N):
    if xf[i] >= 0:
        A =+ yf[i]
    else:
        B =+ yf[i]


omega = 2*np.pi*xf

k = omega/c_0

x_linspace = np.linspace(-100, 100, 100000)

u_x_t = np.zeros(np.shape(x_linspace)[0])

def u(x,t):
    u = 0

    i = int(N/2)

    for n in range(i):
        u =+ np.abs(B*np.exp(-1j*(np.abs(k[n+5])*x + np.abs(omega[n+5])*t)))

    return u

print(omega[5:10])
print(u(0,0))
print(u(0,0.0001))
print(u(0,0.0002))
print(u(0,0.0003))
print(u(0,0.0004))

for i in range(np.shape(x_linspace)[0]):
    u_x_t[i] = u(x_linspace[i],0)

plt.figure('Wave propagation')
plt.plot(x_linspace, u_x_t)
plt.show()