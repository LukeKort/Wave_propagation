# Cálculo da transformada de Fourier

import numpy as np
import matplotlib.pyplot as plt

rho = 7850  # Bar's density
E = 210e9 # Bar's Young Modulos
c_0 = np.sqrt(E/rho) # Bar's propagation velocity

gamma = 2*E

N = 10 # Number of discratization points

A = 1 # Maxium amplitude

bar_len = 5 # Bar lengtg

M = 500 # Mass at the end of the bar

Area = np.pi*(0.0127**2)

a = -1 # Sign limits
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

recoverd_sign = ifft(yf) # Recuperates the original sign from the FFT frequencies

k = 2*np.pi*xf # Wave number

omega = c_0*k # Wave circular frequency

def A_hat(x,t): # Generates the original incident sign
    x = x - c_0*t
    if (x >= -1) and (x <= 0):
        return A*(x + 1)
    if (x > 0) and (x <= 1):
        return A*(-x + 1)
    else:
        return 0
    
x_linspace = np.linspace(-1, 12, 600) # Point over the x axis to calculate u(x,t)

u_x_t = np.zeros(np.shape(x_linspace)[0]) # prealocates the u(x,t) vector

def u_i(x,t): # Calculates the u_incident(x,t) from the espectral form at x and t
    u = 0

    Amp = A_hat(x,t)

    for n in range(N):
        u =+ Amp*np.exp(-1j*(k[n]*x - omega[n]*t))

    return u

def u_t(x,t): # Calculates the u_reflexion(x,t) from the espectral form at x and t\

    u = 0

    Amp_B = A_hat(x,t)

    for n in range(1,N):
        phi = np.sqrt(rho*omega[n]**2/(E-gamma*1j*omega[n]))       
        u =+ Amp_B*np.exp(-1j*(phi*x))

    return u

for j in range(12):
    for i in range(np.shape(x_linspace)[0]): # Calculates the u(x,t) along the x axis at x
        u_x_t[i] = u_t(x_linspace[i], j/c_0) # m/c_0 gives the time need for the wave displace m meters

    # Plots the wave recovered

    title = str('Wave propagation - T = %.4f' %(j/c_0))
    plt.figure(title)
    plt.plot(x_linspace, u_x_t)
    plt.axhline(y = 1, label = 'Original amplitude', linestyle = 'dashed', color = 'g')
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.ylim(0,1.5)
    plt.legend()
    plt.show()