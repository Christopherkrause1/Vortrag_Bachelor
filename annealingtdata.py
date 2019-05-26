import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

N_C0 = 1.3*10**11 #1/cm**3
E_y = 1.33*1.6*10**(-19)    #resulting activation Energy
k_0y = 1.5 * 10**(15)   #frequency factor
g_c = 1.49 * 10**(-2)  #cm**(-1)    Acceptor introduction Rate
g_a = 1.59 * 10**(-2) #cm**(-1)   introduction rate
g_y = 5.16*10**(-2)   #cm**(-1)
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante
E_aa = 1.09 * 1.6* 10**(-19) #j   activation Energy
k_0a = 2.4 *10**(13) #1/s   frequency factor




def N_Y_inf(phi):                                      #longterm annealing amplitude
    return g_y * phi

def tau_Y(T):                                        #Time constant
    return 1/(k_0y *np.exp(-E_y/(k_B*(T+273.15))))

def gett_Y(t, T_s, T):
    timediff_Y = np.zeros(len(t))
    timediff_Y = np.ediff1d(t, to_begin=0)
    T_s = np.roll(T_s, shift=1) # shifting array by one to the right
    T_n = T
    timediff_Y /= tau_Y((T_s+ T_n)/2)
    t_Y = np.zeros(len(t))
    for i in range(0, len(t)):
        t_Y[i] = np.sum(timediff_Y[0:i+1])
    return t_Y


def tau_A(T):                                          #Time constant
    return 1/(k_0a *np.exp(-E_aa/(k_B*(T+273.15))))

def gett_A(t, tau_A0, T):                              #sum of time differences divided by tau(T)
    timediff_A = np.zeros(len(t))
    timediff_A = np.ediff1d(t, to_begin=0)
    tau_A0 = np.roll(tau_A0, shift=1) # shifting array by one to the right, t[1]-t[0] soll ja durch tau[0] geteilt werden
    tau_A1 = tau_A(T)
    timediff_A /= (tau_A0 + tau_A1)/2
    t_A = np.zeros(len(t))
    for i in range(0, len(t)):
        t_A[i] = np.sum(timediff_A[0:i+1])
    return t_A


def N_C(phi):                                          #stable damage
    return N_C0 *(1 - np.exp(-phi)) + g_c * phi

def N_A(t, phi, T):                                    #shortterm annealing
    tau_A0 = tau_A(T)                         #tau_A0 = Array [egal, tau_A(T[0]), tau_A(T[1]),...]
    t_A = gett_A(t, tau_A0, T)                         #Vektor t_1 - t_0/tau_A(0)
    return phi * g_a * np.exp(-t_A)


def N_Y(t, phi, T):                                    #longterm annealing
    T_s = T                             #tau_Y0 = Array [egal, tau_Y(T[0]), tau_Y(T[1]),...]
    t_Y = gett_Y(t, T_s, T)                         #Vektor t_1 - t_0/tau_Y(0)
    return N_Y_inf(phi) * (1- 1/(1 + t_Y))


def N_eff(t, phi, T):                                #Änderung der Dotierungskonzentration
    return N_C(phi) + N_A(t, phi, T) + N_Y(t, phi, T)


#Änderung der effektiven Dotierungskonzentration für R1
t, T_2 = np.genfromtxt('Daten/tdata_1.txt', unpack=True)   #R1 daten


fig, ax1 = plt.subplots()
plt.semilogx(t/60 , T_2, 'r.', label='Temperature', Markersize=6)
#ax1.bar()
#ax1.scatter()
ax1.set_ylabel(r"Temperature / $^{\circ}$C", color = 'red')
ax1.tick_params('y',colors='red')
ax1.set_xlabel("Time / min")
ax1.legend(loc=6)


ax2 = ax1.twinx()
plt.semilogx(t/60, N_eff(t, 5*10**(15), T_2), 'b.', label=r'$\Delta N_{\mathrm{eff}}$ of R1', Markersize=6)
#plt.semilogx(t/60, N_eff(t, 5*10**(15), 80), 'k--', label=r'$\Delta N_{\mathrm{eff}}$ für 80°C', Markersize=6)
ax2.set_ylabel(r"$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $",color='blue')
ax2.tick_params('y',colors='blue')
#ax2.set_yscale('log')
#ax2.scatter()
ax1.grid()
ax2.legend(loc='best')

#plt.gcf().subplots_adjust(bottom=0.18)
#plt.semilogx(t/60, N_eff(t, 5*10**(15), T_2), 'r.', label='Änderung N_eff R1', Markersize=6)
#plt.semilogx(t/60, N_eff(t, 5*10**(15), 80), 'b.', label='Änderung N_eff 80°C', Markersize=6)
#plt.semilogx(t/60, N_C(5*10**(15))+N_A(t, 5*10**(15), T_2) + N_Y(t, 5*10**(15), T_2), 'k-', label='Änderung N_eff R1', Markersize=6)
##plt.semilogx(t/60, N_C(5*10**(15))+N_A(t, 5*10**(15), T_2), 'b-', label='Änderung N_A', Markersize=6)
##plt.semilogx(t/60, N_C(5*10**(15))+N_Y(t, 5*10**(15), T_2), 'g-', label='Änderung N_A', Markersize=6)
#plt.title('Annealingeffekt für R1')
#plt.legend()
#plt.grid()
#plt.xlabel(r't / $\mathrm{min}$')
#plt.ylabel(r'$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $')
plt.savefig('images/annealingtdata.pdf')
plt.clf()




#Änderung der effektiven Dotierungskonzentration für Diode mit Unix Zeiten
#t_unix, T_3 = np.genfromtxt('2018-09-22_11_21_40_Annealingtest_1950.txt', usecols=(0, 2), unpack=True)  #unix daten
#t_s = t_unix-t_unix[0]          #t_s = vergangene Zeit in Sekunden
#
#plt.gcf().subplots_adjust(bottom=0.18)
#plt.semilogx(t_s/60, N_eff(t_s, 1*10**(15), T_3), 'r.', label='Änderung N_eff', Markersize=6)
#plt.semilogx(t_s/60, N_eff(t_s, 1*10**(15), 60), 'b.', label='Änderung N_eff 60°C', Markersize=6)
#plt.title('Annealingeffekt')
#plt.legend()
#plt.grid()
#plt.xlabel(r't / $\mathrm{min}$')
#plt.ylabel(r'$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $')
#plt.savefig('build/annealingunix.pdf')
#plt.clf()



#Zweiter Datensatz mit unix zeiten
#t_unix2, T_4 = np.genfromtxt('2018-09-23_07_40_48_Annealingtest_1950.txt', usecols=(0, 1), unpack=True)  #unix daten
#t_s2 = t_unix2-t_unix2[0]
#
#plt.gcf().subplots_adjust(bottom=0.18)
#plt.semilogx(t_s2/60, N_eff(t_s2, 1*10**(15), T_4), 'r.', label='Änderung N_eff', Markersize=6)
#plt.semilogx(t_s2/60, N_eff(t_s2, 1*10**(15), 60), 'b.', label='Änderung N_eff 60°C', Markersize=6)
#plt.title('Annealingeffekt')
#plt.legend()
#plt.grid()
#plt.xlabel(r't / $\mathrm{min}$')
#plt.ylabel(r'$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $')
#plt.savefig('build/annealingunix_2.pdf')
#plt.clf()
