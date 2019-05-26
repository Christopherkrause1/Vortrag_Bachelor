import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

N_C0 = 1.3*10**11 #1/cm**3
E_y = 1.33*1.6*10**(-19)    #resulting activation Energy
k_0y = 1.5 * 10**(15)*60   #frequnecy factor
g_c = 1.49 * 10**(-2)  #cm**(-1)    Acceptor introduction Rate
g_a = 1.59 * 10**(-2) #cm**(-1)   introduction rate
g_y = 5.16*10**(-2)   #cm**(-1)
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante
E_aa = 1.09 * 1.6* 10**(-19) #j   activation Energy
k_0a = 2.4 *10**(13)*60 #1/min   frequnecy factor


t, phi, T, T_2, T_3, T_4 = np.genfromtxt('Daten/daten.txt', unpack=True)
#t_2, T_5 = np.genfromtxt('Daten/tdata_1.txt', unpack=True)
#t_unix, T_6 = np.genfromtxt('2018-09-22_11_21_40_Annealingtest_1950.txt', usecols=(0, 2), unpack=True)  #unix daten
#t_s = t_unix-t_unix[0]
#Änderung der effektiven Dotierungskonzentration für WE-25k$\Omega$cm


def N_Y_inf(phi):                                      #longterm annealing amplitude
    return g_y * phi

def tau_Y(T):                                          #Time constant
    return 1/(k_0y *np.exp(-E_y/(k_B*T)))

def tau_A(T):                                          #Time constant
    return 1/(k_0a *np.exp(-E_aa/(k_B*T)))

def N_C(phi):                                          #stable damage
    return N_C0 *(1 - np.exp(-phi)) + g_c * phi

def N_A(t, phi, T):                                    #shortterm annealing
    return phi * g_a * np.exp(-t/tau_A(T))

def N_Y(t, phi, T):                                    #longterm annealing
    return N_Y_inf(phi) * (1- 1/(1 + t/tau_Y(T)))


def N_eff(t, phi, T):                                  #Änderung der Dotierungskonzentration
    return N_C(phi) + N_A(t, phi, T) + N_Y(t, phi, T)





#fig, ax1 = plt.subplots()
#plt.semilogx(t_2/60 , T_5, 'r.', label='Temperature', Markersize=6)
#ax1.set_ylabel(r"Temperature / $^{\circ}$C", color = 'red')
#ax1.tick_params('y',colors='red')
#ax1.set_xlabel("Time / min")
#ax1.legend(loc='upper left')
#
#
#ax2 = ax1.twinx()
#plt.semilogx(t_2/60, N_eff(t_2/60, 5*10**(15), T_5+273.15), 'b.', label=r'$\Delta N_{\mathrm{eff}}$ of R1', Markersize=6)
##plt.semilogx(t_2/60, N_eff(t_2/60, 5*10**(15), 80+273.15), 'k--', label=r'$\Delta N_{\mathrm{eff}}$@80°C', Markersize=6)
#ax2.set_ylabel(r"$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $",color='blue')
#ax2.tick_params('y',colors='blue')
#
#ax1.grid()
#ax2.legend(loc='best')
#plt.savefig('build/annealing.pdf')
#plt.clf()

plt.gcf().subplots_adjust(bottom=0.18)
plt.semilogx(t, N_eff(t, 5*10**(15), 60+273.15), 'r.', label=r'$\Delta N_eff@60°C$', Markersize=6)
plt.semilogx(t, N_eff(t, 5*10**(15), 80+273.15), 'b.', label=r'$\Delta N_eff@80°C$', Markersize=6)
plt.legend()
plt.grid()

plt.xlabel(r'Time / $\mathrm{min}$')
plt.ylabel(r'$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $')
plt.savefig('images/annealing.pdf')
plt.clf()
