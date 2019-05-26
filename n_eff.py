from configuration import *

t_1 -= t_1[0]                            #converts unix time stamps to seconds
k_B = 1.38064852 * 10**(-23)                           #Boltzmann constant in J/K



def N_Y_inf(phi):                                      #longterm annealing amplitude
    return g_y * phi

def tau_Y(T):                                          #Time constant of the longterm annealing
    return 1/(k_0y *np.exp(-E_y/(k_B*(T+273.15))))

def gett_Y(t, tau_Y0, T):                    #creating an approximation for t/tau_Y
    timediff_Y = np.zeros(len(t))            #creates an array of zeros to work with
    timediff_Y = np.ediff1d(t, to_begin=0)   #creates array = [0, t[1]-t[0], t[2]-t[1], ...]
    tau_Y0 = np.roll(tau_Y0, shift=1)        #shifting tau_Y array by one to the right
    tau_Y1 = tau_Y(T)
    timediff_Y /= (tau_Y0+ tau_Y1)/2         #dividing each element by the mean of 2 adjacent tau_Y elements
    t_Y = np.zeros(len(t))                   #create an array of zeros to work with
    for i in range(0, len(t)):
        t_Y[i] = np.sum(timediff_Y[0:i+1])   #writes in each element the sum of the cooresponding timediff_Y values
    return t_Y                               #now looks like [0, 0 + t[1]-t[0]/(tau_Y[0] + tau_Y[1])/2, ...]


def tau_A(T):                                #Time constant of the short term annealing
    return 1/(k_0a *np.exp(-E_aa/(k_B*(T+273.15))))

def gett_A(t, tau_A0, T):                   #creating an approximation for t/tau_A
    timediff_A = np.zeros(len(t))           #creates an array of zeros to work with
    timediff_A = np.ediff1d(t, to_begin=0)  #creates array = [0, t[1]-t[0], t[2]-t[1], ...]
    tau_A0 = np.roll(tau_A0, shift=1)       #shifting tau_A array by one to the right
    tau_A1 = tau_A(T)
    timediff_A /= (tau_A0 + tau_A1)/2       #dividing each element by the mean of 2 adjacent tau_A elements
    t_A = np.zeros(len(t))                  #create an array of zeros to work with
    for i in range(0, len(t)):
        t_A[i] = np.sum(timediff_A[0:i+1])  #writes in each element the sum of the cooresponding timediff_Y values
    return t_A                              #now looks like [0, 0 + t[1]-t[0]/(tau_A[0] + tau_A[1])/2, ...]


def N_C(phi):                                        #stable damage
    return N_C0 *(1 - np.exp(-c * phi)) + g_c * phi

def N_A(t, phi, T):                                  #shortterm annealing
    tau_A0 = tau_A(T)                                #assigning array for the approximation
    t_A = gett_A(t, tau_A0, T)                       #assigning name to the new approximated time
    return phi * g_a * np.exp(-t_A)


def N_Y(t, phi, T):                                  #longterm annealing
    tau_Y0 = tau_Y(T)                                #assigning array for the approximation
    t_Y = gett_Y(t, tau_Y0, T)                       #assigning name to the new approximated time
    return N_Y_inf(phi) * (1- 1/(1 + t_Y))


def interpolation_t(t, T):                                            #linear interpolation of the time for more data
    t_int = np.array(t[0])                                            #create new time starting with arrays of zeros
    T_min = min(T)                                                    #number of intervalls depend on minimal temperature
    for i in range(1, len(T)):
        n = math.ceil((x_int*abs(T[i-1]- T_min) + y_int))             #function for the number of intervalls
        for j in range(1, n+1):
            t_int = np.append(t_int, t[i-1] + abs(t[i-1]-t[i])/n *j)  #new interpolated times (includes initial times)
    return t_int;


def interpolation_T(t, T):                                            #linear interpolation of the temperature for more data
    T_int = np.array(T[0])                                            #create new temperature starting with arrays of zeros
    T_min = min(T)                                                    #number of intervalls depend on minimal temperature
    for i in range(1, len(T)):
        n = math.ceil((x_int*abs(T[i-1]- T_min) + y_int))             #function for the number of intervalls
        for j in range(1, n+1):
            T_int = np.append(T_int, T[i-1] + (T[i]-T[i-1])/(n) * (j))#new interpolated temperatures (includes initial ones)
    return T_int




def N_eff(t, phi, T):  #change of the doping concentration
    return N_C(phi) + N_A(t, phi, T) + N_Y(t, phi, T)




#function that plots n_eff
def plot_N_eff(t, phi, T):
    fig, ax1 = plt.subplots()
    plt.semilogx(interpolation_t(t, T)/60 , interpolation_T(t, T), 'g.', label='interpolated temperature', Markersize=6)
    plt.semilogx(t/60 , T, 'r.', label='measured temperature', Markersize=6)
    ax1.set_ylabel(r"Temperature / $^{\circ}$C", color = 'red')
    ax1.tick_params('y',colors='red')
    ax1.set_xlabel("Time / min")
    ax1.legend(loc=6)


    ax2 = ax1.twinx()
    plt.semilogx(interpolation_t(t, T)/60, N_eff(interpolation_t(t, T), phi, interpolation_T(t, T)), 'b.', label=r'interpolated $\Delta N_{\mathrm{eff}}$', Markersize=6)
    plt.semilogx(t/60, N_eff(t, phi, T), 'k.', label=r'$\Delta N_{\mathrm{eff}}$')
    ax2.set_ylabel(r"$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $",color='blue')
    ax2.tick_params('y',colors='blue')
    ax1.grid()
    ax2.legend(loc='best')
    plt.savefig('images/interpolationtdata.pdf')
    plt.clf()

plot_N_eff(t_1, phi, T_1)    #function of the doping concentration

#zweiter Datensatz
#fig, ax1 = plt.subplots()
#plt.semilogx(interpolation_t(t_3, T_3)/60 , interpolation_T(t_3, T_3), 'r.', label='interpolierte Temperatur', Markersize=6)
#plt.semilogx(t_3/60 , T_3, 'g.', label='Temperatur', Markersize=6)
#ax1.set_ylabel(r"Temperatur / $^{\circ}$C", color = 'red')
#ax1.tick_params('y',colors='red')
#ax1.set_xlabel("Zeit / min")
#ax1.legend(loc=6)
#
#
#ax2 = ax1.twinx()
#plt.semilogx(interpolation_t(t_3, T_3)/60, N_eff(interpolation_t(t_3, T_3), phi, interpolation_T(t_3, T_3)), 'b.', label=r'$\Delta N_{\mathrm{eff}}$ mit interpolation', Markersize=6)
#plt.semilogx(t_3/60, N_eff(t_3, 5*10**(15), T_3), 'k.', label=r'$\Delta N_{\mathrm{eff}}$ ohne interpolation', Markersize=6)
#ax2.set_ylabel(r"$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $",color='blue')
#ax2.tick_params('y',colors='blue')
#ax1.grid()
#ax2.legend(loc='best')
#plt.savefig('build/interpolationmareike2.pdf')
#plt.clf()
#
#
##erster Datensatz R3
#fig, ax1 = plt.subplots()
#plt.semilogx(interpolation_t(t_4, T_4)/60 , interpolation_T(t_4, T_4), 'r.', label='interpolierte Temperatur', Markersize=6)
#plt.semilogx(t_4/60 , T_4, 'g.', label='Temperatur', Markersize=6)
#ax1.set_ylabel(r"Temperatur / $^{\circ}$C", color = 'red')
#ax1.tick_params('y',colors='red')
#ax1.set_xlabel("Zeit / min")
#ax1.legend(loc=6)
#
#
#ax2 = ax1.twinx()
#plt.semilogx(interpolation_t(t_4, T_4)/60, N_eff(interpolation_t(t_4, T_4), phi, interpolation_T(t_4, T_4)), 'b.', label=r'$\Delta N_{\mathrm{eff}}$ mit interpolation', Markersize=6)
#plt.semilogx(t_4/60, N_eff(t_4, 5*10**(15), T_4), 'k.', label=r'$\Delta N_{\mathrm{eff}}$ ohne interpolation', Markersize=6)
#ax2.set_ylabel(r"$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $",color='blue')
#ax2.tick_params('y',colors='blue')
#ax1.grid()
#ax2.legend(loc='best')
#plt.savefig('build/interpolationmareikeR3_1.pdf')
#plt.clf()
