"""
Purpose:
"""

def euler_integration(func, xvt_init, dt, n):
    x = np.zeros(n); x[0] = xvt_init[0]
    v = np.zeros(n); v[0] = xvt_init[1]
    t = np.zeros(n); t[0] = xvt_init[2]

    for i in range(n-1):
        v[i+1] = v[i] + func(x[i],v[i],t[i])*dt 
        x[i+1] = x[i] + v[i]*dt
        t[i+1] = t[i] + dt

    return x,v,t

def simple_predcorr_integration(func, xvt_init, dt, n):
    x = np.zeros(n); x[0] = xvt_init[0] #position-array
    v = np.zeros(n); v[0] = xvt_init[1] #velocity-array
    t = np.zeros(n); t[0] = xvt_init[2] #time-array
    dt_2 = dt/2.0 #pre-calculation
    
    for i in range(n-1):
        t[i+1] = t[i] + dt

        #make prediction with euler
        a_i = func(x[i],v[i],t[i])
        v[i+1] = v[i] + a_i*dt 
        x[i+1] = x[i] + v[i]*dt

        #make correction to prediction with the trapezoid rule
        a_ii = func(x[i+1],v[i+1],t[i+1])
        v[i+1] = v[i] + (a_i+a_ii)*dt_2
        x[i+1] = x[i] + (v[i]+v[i+1])*dt_2

    return x,v,t

def euler_cromer_integration(func, xvt_init, dt, n):
    x = np.zeros(n); x[0] = xvt_init[0]
    v = np.zeros(n); v[0] = xvt_init[1]
    t = np.zeros(n); t[0] = xvt_init[2]

    for i in range(n-1):
        v[i+1] = v[i] + func(x[i],v[i],t[i])*dt 
        x[i+1] = x[i] + v[i+1]*dt
        t[i+1] = t[i] + dt

    return x,v,t

def midpoint_integration(func, xvt_init, dt, n):
    x = np.zeros(n); x[0] = xvt_init[0]
    v = np.zeros(n); v[0] = xvt_init[1]
    t = np.zeros(n); t[0] = xvt_init[2]

    for i in range(n-1):
        x_half = x[i] + v[i]*dt/2.0
        v_half = v[i] + func(x[i],v[i],t[i])*dt/2.0
        a_half = func(x_half, v_half, t[i] + dt/2.0)

        x[i+1] = x[i] + v_half*dt
        v[i+1] = v[i] + a_half*dt 
        
        t[i+1] = t[i] + dt

    return x,v,t

def rungekutta4_integration(func, xvt_init, dt, n):
    x = np.zeros(n); x[0] = xvt_init[0]
    v = np.zeros(n); v[0] = xvt_init[1]
    t = np.zeros(n); t[0] = xvt_init[2]

    for i in range(n-1):
        dt_2 = dt/2.0
        dt_6 = dt/6.0
        
        x1 = x[i]
        v1 = v[i]
        a1 = func(x1,v1,t[i])
        
        x2 = x[i] + v1*dt_2
        v2 = v[i] + a1*dt_2
        a2 = func(x2,v2,t[i])

        x3 = x[i] + v2*dt_2
        v3 = v[i] + a2*dt_2
        a3 = func(x3,v3,t[i])

        x4 = x[i] + v3*dt
        v4 = v[i] + a3*dt
        a4 = func(x4,v4,t[i])
        
        v[i+1] = v[i] + (v1 + 2.0*v2 + 2.0*v3 + v4)*dt_6
        x[i+1] = x[i] + (x1 + 2.0*x2 + 2.0*x3 + x4)*dt_6
        t[i+1] = t[i] + dt

    return x,v,t



def integrate_time(func, init, timearray, method):
    arraylength = len(timearray)
    num_integration = len(init)
    if method.lower() == "euler":
        timestep = timearray[1] - timearray[0]
        init = list(init) + [timearray[0]] #add initial time to list of init
        print "Warning!"
        print "method 'euler_twostep_integration()' assumes"
        print "linearly spaced time"
        x,v,t = euler_integration(func, xvt_init=init,
                                  dt=timestep, n=arraylength)
    elif method.lower() == "prediction/correction" or method.lower() == "predcorr":
        timestep = timearray[1] - timearray[0]
        init = list(init) + [timearray[0]] #add initial time to list of init
        x,v,t = simple_predcorr_integration(func, xvt_init=init,
                                            dt=timestep, n=arraylength)
    elif method.lower() == "eulercromer":
        timestep = timearray[1] - timearray[0]
        init = list(init) + [timearray[0]] #add initial time to list of init
        x,v,t = euler_cromer_integration(func, xvt_init=init,
                                  dt=timestep, n=arraylength)
    elif method.lower() == "midpoint":
        timestep = timearray[1] - timearray[0]
        init = list(init) + [timearray[0]] #add initial time to list of init
        x,v,t = midpoint_integration(func, xvt_init=init,
                                  dt=timestep, n=arraylength)
    elif method.lower() == "rungekutta4" or method.lower() == "rk4":
        timestep = timearray[1] - timearray[0]
        init = list(init) + [timearray[0]] #add initial time to list of init
        x,v,t = rungekutta4_integration(func, xvt_init=init,
                                  dt=timestep, n=arraylength)
    return x,v,t
        
def test_scenario(timestep_per_pi, int_method):
    """
    Test integration-methods for a given timestep.
    Consider a mass on a spring, with the oscillatory force F=-kx=ma.
    Solve the 1D system of ODEs; dv/dt = a; dx/dt = v
    k = m = v0 = 1; x0 = 0
    Analytical solution: x(t) = x_m sin(wt); v(t) = x_m w cos(wt); x_m = w = 1
    Total energy: E = 1/2 mv^2 + 1/2 kx^2
    
    The analysis will be executed with 2 diff. methods (Euler and Pred-Corr),
    for varying degrees of timestep
    """

    #determine BC and IC
    x0 = 0.0 #init pos
    v0 = 1.0 #init vel
    t0 = 0.0 #start-time
    tn = 4.0*np.pi #end-time
    tau = timestep_per_pi*np.pi #timesteps
    n = (tn-t0)/tau + 1 #number of timesteps
    
    time = np.linspace(t0, tn, n) #time-array

    #acceleration of point particle with k=m=1
    acc1 = lambda x,v,t: -1.0*x #function must take three arguments!

    pos, vel, time = integrate_time(func=acc1,
                                    init=(x0,v0),
                                    timearray=time,
                                    method=int_method)

    #analytical solutions
    pos_an = np.sin(time)
    vel_an = np.cos(time)

    return time, pos, pos_an, vel, vel_an


def total_energy(x,v):
    Ep = 0.5*x*x #k=1
    Ek = 0.5*v*v #m=1
    return Ep + Ek

def relative_error(e, e0):
    return np.abs((e-e0)/e0) #relative error in energy


if __name__ =='__main__':
    """
    Purpose:
    """
    import matplotlib.pyplot as pl

    epsilon = {}
    
    for method in ["euler", "eulercromer", "rungekutta4"]:
        timesteps = 10**np.linspace(-1,-5,5)
        epsilon[method] = np.zeros((3,2))
        for dt in timesteps:
            t, x_num, x_an, v_num, v_an = test_scenario(timestep_per_pi=dt, int_method=method)
            E_tot_start = total_energy(x_num[0],v_num[0])
            E_tot_end = total_energy(x_num[-1],v_num[-1])
            epsilon_current = relative_error(E_tot_end, E_tot_start)
            epsilon[method][0,:] = np.array([epsilon_current, dt])
            

    pl.figure(); pl.grid(True)
    for key in epsilon.keys():
        e = epsilon[key][:,0] #all errors for a given method
        dt = epsilon[key][:,1] #timesteps of errors
        pl.plot(dt,e, '.-', label=key)

    pl.xlabel("timesteps " + r"[$\pi$]")
    pl.ylabel("relative error " + r"$\left|\frac{\epsilon_n-\epsilon_0}{\epsilon_0}\right|$")

    pl.xscale("log")
    pl.yscale("log")
    pl.legend(loc='best')
    #pl.savefig()
    pl.show()
    
