from configuration import *
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante




def N_Y_inf(phi):                                      #longterm annealing amplitude
    return g_y * phi

def tau_Y(T):                                        #Time constant
    return 1/(k_0y *np.exp(-E_y/(k_B*(T+273.15))))

def gett_Y(t, tau_Y0, T):
    timediff_Y = np.zeros(len(t))
    timediff_Y = np.ediff1d(t, to_begin=0)
    tau_Y0 = np.roll(tau_Y0, shift=1) # shifting array by one to the right
    tau_Y1 = tau_Y(T)
    timediff_Y /= (tau_Y0+ tau_Y1)/2
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
    t_A = gett_A(t, tau_A0, T)                   #Vektor t_1 - t_0/tau_A(0)
    return phi * g_a * np.exp(-t_A)


def N_Y(t, phi, T):                                    #longterm annealing
    tau_Y0 = tau_Y(T)                             #tau_Y0 = Array [egal, tau_Y(T[0]), tau_Y(T[1]),...]
    t_Y = gett_Y(t, tau_Y0, T)                         #Vektor t_1 - t_0/tau_Y(0)
    return N_Y_inf(phi) * (1- 1/(1 + t_Y))


def N_eff(t, phi, T):                                #Änderung der Dotierungskonzentration
    return N_C(phi) + N_A(t, phi, T) + N_Y(t, phi, T)




new_t = np.array(t_1[0])
new_T = np.array(T_1[0])

T_min = min(T_1)

for i in range(1, len(T_1)):
    n = math.ceil((x_int*abs(T_1[i-1]- T_min) + y_int))
    for j in range(1, n+1):
        new_T = np.append(new_T, T_1[i-1] + (T_1[i]-T_1[i-1])/(n) * (j))
        new_t = np.append(new_t, t_1[i-1] + abs(t_1[i-1]-t_1[i])/n *j)



#fig, ax1 = plt.subplots()
#plt.semilogx(new_t/60 , new_T, 'r.', label='interpolated Temperature', Markersize=6)
#plt.semilogx(t_1/60 , T_1, 'g.', label='Temperature', Markersize=6)
#ax1.set_ylabel(r"Temperatur / $^{\circ}$C", color = 'red')
#ax1.tick_params('y',colors='red')
#ax1.set_xlabel("Zeit / min")
#ax1.legend(loc=6)
#
#
#ax2 = ax1.twinx()
#plt.semilogx(new_t/60, N_eff(new_t, 5*10**(15), new_T), 'b.', label=r'$\Delta N_{\mathrm{eff}}$ with interpolation', Markersize=6)
#plt.semilogx(t_1/60, N_eff(t_1, 5*10**(15), T_1), 'k.', label=r'$\Delta N_{\mathrm{eff}}$ without interpolation', Markersize=6)
##plt.semilogx(new_t/60, N_C(5*10**(15))+N_Y(new_t, 5*10**(15), new_T), 'y-', label=r'$\Delta N_{\mathrm{eff}}$ ohne interpolation', Markersize=6)
##plt.semilogx(new_x/60, N_eff(new_x, 5*10**(15), 80), 'k--', label=r'$\Delta N_{\mathrm{eff}}$ für 80°C', Markersize=6)
#ax2.set_ylabel(r"$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $",color='blue')
#ax2.tick_params('y',colors='blue')
#ax1.grid()
#ax2.legend(loc='best')
##plt.xlim(0, 2*10**2)
##plt.savefig('images/interpolationtdata.pdf')
#plt.clf()

plt.gcf().subplots_adjust(bottom=0.18)
plt.semilogx(new_t/60, new_T, 'g.', label='interpolated Temperature', Markersize=6)
plt.semilogx(t_1/60, T_1, 'r.', label='measured Temperature', Markersize=6)
#plt.semilogx(new_t/60, N_C(5*10**(15))+N_Y(new_t, 5*10**(15), new_T), 'y-', label=r'$\Delta N_{\mathrm{eff}}$ ohne interpolation', Markersize=6)
#plt.semilogx(new_x/60, N_eff(new_x, 5*10**(15), 80), 'k--', label=r'$\Delta N_{\mathrm{eff}}$ für 80°C', Markersize=6)
plt.ylabel(r"$T/ °C$")
plt.xlabel(r"Time$/ min$")
plt.grid()
plt.legend(loc='best')
#plt.xlim(0, 2*10**2)
plt.savefig('images/interpolationtemperatur.pdf')
plt.clf()
