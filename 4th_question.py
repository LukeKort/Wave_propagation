# Cálculo da transformada de Fourier

import numpy as np
import matplotlib.pyplot as plt

rho = 7850  # Bar's density
E = 210e9 # Bar's Young Modulos
c_0 = np.sqrt(E/rho) # Bar's propagation velocity

N = 1000 # Number of discratization points

A = 1 # Maxium amplitude

a = -2 # Sign limits
b = 0
c = 1

x_disc = np.linspace(-2,1,N) # Discretizes the sign along the y axis

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

# Recovery of the sign

recoverd_sign = ifft(yf) # Recuperates the original sign from the FFT frequencies

# Calculate the wave propagation

k = 2*np.pi*xf # Wave number

omega = c_0*k # Wave circular frequency

def A_hat(x,t): # Generates the original sign
    x = x + c_0*t
    if (x >= -2) and (x <= 0):
        return A*(x/2 + 1)
    if (x > 0) and (x <= 1):
        return A*(-x + 1)
    else:
        return 0

x_linspace = np.linspace(-6, 6, 600) # Point over the x axis to calculate u(x,t)
 
u_x_t = np.zeros(np.shape(x_linspace)[0]) # prealocates the u(x,t) vector

def u(x,t):  # Calculates the u(x,t) from the espectral form at x and t
    
    u = 0

    Amp_A = A_hat(x,-t)


    for n in range(N):
        u =+ Amp_A*np.exp(-1j*(k[n]*x - omega[n]*t))

    return u

def u_love(x,t):  # Calculates the u(x,t) from the espectral form at x and t
    
    u = 0

    gamma = 1

    Amp_A = A_hat(x,-t)
    
    
    for i in range(4):
        u =+ Amp_A*np.exp(1j*gamma*(x - c_0*t))

    return u

for j in range(10):
    for i in range(np.shape(x_linspace)[0]): # Calculates the u(x,t) along the x axis at x
        u_x_t[i] = u_love(x_linspace[i],j/c_0) # m/c_0 gives the time needed for the wave displace m meters

    # Plots the wave recovered

    plt.figure('Wave propagation')
    plt.plot(x_linspace, u_x_t)
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.show()