#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def process_sounding(fname):
    p,t,td,rh,hgt = np.loadtxt(fname, skiprows=15, usecols=(1,2,3,4,14), unpack=True)
    
    # mask levels with missing p, t, td, rh, or hgt data
    coldpool_mask = (p < 9000) & (t < 900) & (td < 900) & (rh < 9000) & (hgt < 90000)
    p,t,td,rh,hgt = p[coldpool_mask], t[coldpool_mask], td[coldpool_mask], rh[coldpool_mask], hgt[coldpool_mask]

    # compute virtual temperature
    es = np.power(10, (9.4041-(2354/(t+273.15))))
    e = (rh/100.0)*es
    w = 0.622*(e/(p-e))
    tv = (t+273.15)*(1+0.61*w)

    # compute density
    rho = (p*100)/(287*tv)

    return (p, t, tv, hgt, rho)

# process each sounding using handy function above
p, t, tv, hgt, rho = process_sounding('sounding_coldpool.txt')
p_env, t_env, tv_env, hgt_env, rho_env = process_sounding('sounding_environment.txt')

# set some variables
b_all = []
integrand_1, integrand_2 = 0,0

# loop through each level in cold pool sounding
for i in range(len(p)):
    # find environmental level closest to this cold pool level
    j = np.argmin(np.abs(hgt_env - hgt[i]))

    # compute buoyancy at this level, using closest level in environmental sounding
    b = ((tv[i]-tv_env[j])/tv_env[j])*9.81
    b_all.append(b)

    # integrate levels up to level of neutral buoyancy
    # use centered difference to estimate area in each layer
    if b < 0 and i>0:
        # cold pool assuming boussinesq
        integrand_1 += 0.5*(b_all[i]+b_all[i-1])*(hgt[i]-hgt[i-1])

        # cold pool assuming anelatic
        integrand_2 += 0.5*(b_all[i]+b_all[i-1])*rho_env[j]*(hgt[i]-hgt[i-1])

# print analytic cold pool speed assuming boussinesq atmosphere (M&R eqn 5.45)
print np.sqrt(-2*integrand_1)

# print analytic cold pool speed assuming anelastic atmosphere
print np.sqrt(-(2/rho_env[0])*integrand_2)

# print cold pool speed using surface buoyancy estimate (M&R eqn 5.46)
coldpooldepth = 3600
print np.sqrt(-((tv[0]-tv_env[0])/tv_env[0])*9.81*coldpooldepth)

# plot buoyancy profile
plt.plot(b_all, -p)
plt.xlim((-0.4,0.4))
plt.show()
