import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


#current related damage rate
a_I = 1.23 * 10**(-17) #A/cm
k_0I = 1.2 * 10**(13)*60 #1/min
E_I = 1.11 * 1.6 * 10**(-19) #j
b = 3.07*10**(-18)    #A/cm
t_0 = 1 #min
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante

t, phi, T, T_2, T_3, T_4 = np.genfromtxt('Daten/daten.txt', unpack=True)
t_5, T_5 = np.genfromtxt('Daten/tdata_1.txt', unpack=True)
#t_1, T_1 = np.genfromtxt('Daten/tdata.txt', unpack=True)


def tau_I(T):                                     #time constant
    return 1/(k_0I* np.exp(-E_I/(k_B*T)))

#def gett_I(t, tau_I0, T):
#    timediff_I = np.zeros(len(t))
#    timediff_I = np.ediff1d(t, to_begin=0)
#    tau_I0 = np.roll(tau_I0, shift=1) # shifting array by one to the right
#    tau_I1 = tau_I(T)
#    timediff_I /= (tau_I0 + tau_I1)/2
#    t_I = np.zeros(len(t))
#    for i in range(0, len(t)):
#        t_I[i] = np.sum(timediff_I[0:i+1])
#    return t_I

def a_0(T):                                       #part of the longterm annealing
    return -8.9*10**(-17) + 4.6*10**(-14) * 1/T

def damage(t, T):
    #tau_I0 = tau_I(T)                                  #tau_I0 = Array [egal, tau_I(T[0]), tau_I(T[1]),...]
    #t_I = gett_I(t, tau_I0, T)                            #damage rate
    return (a_I * np.exp(-t/tau_I(T)) + a_0(T) - b * np.log((t+10**(-50))/t_0))




fig, ax1 = plt.subplots()
plt.semilogx(t_5/60 , T_5, 'r.', label='Temperature', Markersize=6)
ax1.set_ylabel(r"Temperature / $^{\circ}$C", color = 'red')
ax1.tick_params('y',colors='red')
ax1.set_xlabel("Time / min")
ax1.legend(loc='upper left')


ax2 = ax1.twinx()
plt.semilogx(t_5/60, damage(t_5/60, T_5+273.15), 'b.', label=r'$\Delta N_{\mathrm{eff}}$ of R1', Markersize=6)
#plt.semilogx(t_2/60, N_eff(t_2/60, 5*10**(15), 80+273.15), 'k--', label=r'$\Delta N_{\mathrm{eff}}$@80°C', Markersize=6)
ax2.set_ylabel(r"$\alpha $ /$\mathrm{A cm^{-1}} $",color='blue')
ax2.tick_params('y',colors='blue')
plt.ylim(0, 1*10**(-16))
ax1.grid()
ax2.legend(loc='best')
plt.xlabel(r'Time / $\mathrm{min}$')
plt.savefig('images/damage.pdf')
plt.clf()
