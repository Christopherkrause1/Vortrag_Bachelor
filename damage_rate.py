from configuration import *
t_1 -= t_1[0]                            #converts unix time stamps to seconds
k_B = 1.38064852 * 10**(-23)             #Boltzmann constant in J/K
t_0 = 60                                 #s



def tau_I(T):                            #time constant
    return 1/(k_0I* np.exp(-E_I/(k_B*T)))

def gett_I(t, tau_I0, T):                   #creates an approximation for t/tau_I
    timediff_I = np.zeros(len(t))           #creates an array of zeros to work with
    timediff_I = np.ediff1d(t, to_begin=0)  #creates array = [0, t[1]-t[0], t[2]-t[1], ...]
    tau_I0 = np.roll(tau_I0, shift=1)       #shifting tau_I array by one to the right
    tau_I1 = tau_I(T)
    timediff_I /= (tau_I0 + tau_I1)/2       #dividing each element by the mean of 2 adjacent tau_I elements
    t_I = np.zeros(len(t))                  #create an array of zeros to work with
    for i in range(0, len(t)):
        t_I[i] = np.sum(timediff_I[0:i+1])  #writes in each element the sum of the cooresponding timediff_I values
    return t_I                              #now looks like [0, 0 + t[1]-t[0]/(tau_I[0] + tau_I[1])/2, ...]


def a_02():                                  #temperature independent
    return a_0 + (b_0/ T_ref)

def theta(T):                                #scaling factor for the time
    return np.exp(-E_I2/k_B *(1/T - 1/T_ref))


def gett_theta(t, theta_0, T):                        #creates an approximation for theta*t
    timediff_theta = np.zeros(len(t))                 #creates an array of zeros to work with
    timediff_theta = np.ediff1d(t, to_begin=0)        #creates array = [0, t[1]-t[0], t[2]-t[1], ...]
    theta_0 = np.roll(theta_0, shift=1)               #shifting theta array by one to the right
    theta_1 = theta(T)
    timediff_theta *= (theta_0 + theta_1)/2           #dividing each element by the mean of 2 adjacent theta elements
    timediff_theta[0]=10**(-90)                       #to avoid dividing by zero in the logarithm
    t_theta = np.zeros(len(t))                        #create an array of zeros to work with
    for i in range(0, len(t)):
        t_theta[i] = np.sum(timediff_theta[0:i+1])    #writes in each element the sum of the cooresponding timediff values
    return t_theta                                    #now looks like [0, 0 + t[1]-t[0]*(theta[0] + theta[1])/2, ...]

def damage(t, T):                                      #damage rate
    tau_I0 = tau_I(T)                                  #assigning array for the t/tau_I approximation
    t_I = gett_I(t, tau_I0, T)                         #assigning name to the new approximated time
    theta_0 = theta(T)                                 #assigning array theta*t for the approximation
    t_theta = gett_theta(t, theta_0, T)                #assigning name for the approximation
    return a_I * np.exp(-t_I) + a_02() - beta * np.log(t_theta /t_0)

def interpolation_t(t, T):                                #linear interpolation of the time for more data
    t_int = np.array(t[0])                                #create new time starting with arrays of zeros
    T_min = min(T)                                        #number of intervalls depend on minimal temperature
    for i in range(1, len(T)):
        n = math.ceil((x_int*abs(T[i-1]- T_min) + y_int)) #function for the number of intervalls
        for j in range(1, n+1):
            t_int = np.append(t_int, t[i-1] + abs(t[i-1]-t[i])/n *j)
    return t_int;                                         #new interpolated times (includes initial times)


def interpolation_T(t, T):                                #linear interpolation of the temperature for more data
    T_int = np.array(T[0])                                #create new temperature starting with arrays of zeros
    T_min = min(T)                                        #number of intervalls depend on minimal temperature
    for i in range(1, len(T)):
        n = math.ceil((x_int*abs(T[i-1]- T_min) + y_int)) #function for the number of intervalls
        for j in range(1, n+1):
            T_int = np.append(T_int, T[i-1] + (T[i]-T[i-1])/(n) * (j))
    return T_int                                          #new interpolated temperatures (includes initial temperatures)


#function that plots the damage rate
def plot_damage_rate(t, T):
    fig, ax1 = plt.subplots()
    plt.semilogx(interpolation_t(t, T)/60 , interpolation_T(t, T), 'r.', label='interpolated temperature', Markersize=6)
    plt.semilogx(t/60 , T, 'g.', label='measured temperature', Markersize=6)
    ax1.set_ylabel(r"Temperature / $^{\circ}$C", color = 'red')
    ax1.tick_params('y',colors='red')
    ax1.set_xlabel("Time / min")
    ax1.legend(loc='best')


    ax2 = ax1.twinx()

    #plt.semilogx(t_1/60 , damage(t_1, 49+273.15), 'b.', label='Schadensrate 49°C', Markersize=6)
    plt.semilogx(interpolation_t(t, T)/60 , damage(interpolation_t(t, T), interpolation_T(t, T)+273.15), 'b.', label='interpolated damage rate', Markersize=6)
    plt.semilogx(t/60 , damage(t, T+273.15), 'k.', label='damage rate', Markersize=6)
    ax2.set_ylabel(r"$\alpha  / \mathrm{A cm^{-1}} $",color='blue')
    ax2.tick_params('y',colors='blue')
    plt.ylim(0.1*10**(-16), 0.9*10**(-16))
    ax1.grid()
    ax2.legend(loc='lower center')
    #plt.show()
    plt.savefig('images/damageinterpolation.pdf')
    plt.clf()

#function of the damage rate
plot_damage_rate(t_1, T_1)
