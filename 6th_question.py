import numpy as np
import matplotlib.pyplot as plt

rho = 7850  # Bar's density
E = 210e9 # Bar's Young Modulos
c_0 = np.sqrt(E/rho) # Bar's propagation velocity

# Fist case - 1 element bar

L = 0.3 # Bar length

r = 0.0125 # Bar's radio

A = np.pi*(r**2) # Area

P = 100 # Load at right end

M = rho*A*L # Bar's mass

I = 0.5*M*r**2 # Moment of inertia

N = 1000

def u_2(w): # Calculates u_2 for given angular frequency

    k = w/c_0 # wave number

    u = P/((E*A/L)*((1j*k*L)/(1-np.exp(-1j*2*k*L)))*(1+np.exp(-1j*2*k*L)))

    return u

u_w = np.zeros(4) # Stores u_2 for given u_w omega
u_w_plus = np.zeros(N)


w = np.zeros(4) # Stores the frequencies
w_plus = np.zeros(len(u_w_plus)) # Stores the frequencies

def w_n(n): # Calculates the n natural frequency

    n = n + 1

    f_1 = 2*np.pi*0.162*(2*r/L**2)*np.sqrt(E/rho)

    if n == 1:
    
        w = f_1

    else:

        w = 2*np.pi*(2.81*(n-0.5)**2)*f_1

    # w = (n**2)*np.sqrt((E*I)/(rho*A*L**4))

    return w

# n = [1001,1002,1003,1004]
n = [1,2,3,4]

n_plus = np.linspace(n[0],n[-1],len(u_w_plus))

for i in range(4):
    w[i] = w_n(n[i])

print(w)

for i in range(4):
    u_w[i] = np.abs(u_2(w[i]))

for i in range(len(u_w_plus)):
    w_plus[i] = w_n(n_plus[i])

for i in range(len(u_w_plus)):
    u_w_plus[i] = np.abs(u_2(w_plus[i]))

title = str('1st to 4th Natural Frequency')
# title = str('1001th to 1004th Natural Frequency')
plt.figure(title)
plt.plot(w_plus,u_w_plus)
plt.scatter(w,u_w, color = 'green')
plt.yscale('log')
plt.xlabel('\u03C9 [rads/s]')
plt.ylabel('u(\u03C9) [m/N]')
plt.title(title)
plt.show()