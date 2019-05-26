import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


#current related damage rate
a_I = 1.23 * 10**(-17) #A/cm
k_0I = 1.2 * 10**(13) #1/s
E_I = 1.11 * 1.6 * 10**(-19) #j
b = 3.07*10**(-18)    #A/cm
t_0 = 1 #min
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante
a_I = 1.23 * 10**(-17)       #amplitude in A/cm
a_0 = -8.9*10**(-17)         #fit parameter in A/cm
k_0I = 1.2 * 10**(13)        #fit parameter 1/s
E_I = 1.11 * 1.6 * 10**(-19) #fit parameter j
E_I2 = 1.3 * 1.6 * 10**(-19) #fit parameter j
beta = 3.07*10**(-18)        #fit parameter A/cm
b_0 = 4.6*10**(-14)          #fit parameter in A*K/cm
T_ref = 353.15               #reference temperature in kelvin

t_5, T_5 = np.genfromtxt('Daten/tdata_1.txt', unpack=True)
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

fig, ax1 = plt.subplots()
plt.semilogx(t_5/60 , T_5, 'r.', label='Temperature', Markersize=6)
ax1.set_ylabel(r"Temperature / $^{\circ}$C", color = 'red')
ax1.tick_params('y',colors='red')
ax1.set_xlabel("Time / min")
ax1.legend(loc='upper left')


ax2 = ax1.twinx()
plt.semilogx(t_5/60, damage(t_5/60, T_5+273.15), 'b.', label=r'$\Delta N_{\mathrm{eff}}$ of R1', Markersize=6)
ax2.set_ylabel(r"$\alpha $ /$\mathrm{A cm^{-1}} $",color='blue')
ax2.tick_params('y',colors='blue')
plt.ylim(0, 1*10**(-16))
ax1.grid()
ax2.legend(loc='best')
plt.xlabel(r'Time / $\mathrm{min}$')
plt.ylim(0.4*10**(-16), 1*10**(-16))
plt.savefig('images/damagekorrektur.pdf')
plt.clf()
