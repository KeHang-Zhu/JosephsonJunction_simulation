#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import kwant
import tkwant



def make_system(length):
    
    def onsite_potential(site, time):
        return 1 + v(time)  # one is the static onsite element

    # system building
    lat = kwant.lattice.square(a=1, norbs=1)
    syst = kwant.Builder()

    # central scattering region
    syst[(lat(x, 0) for x in range(length))] = 1#dont understand
    syst[lat.neighbors()] = -1#dont understand
    # time dependent onsite-potential V(t) at leftmost site
    syst[lat(0, 0)] = onsite_potential

    # add leads
    sym = kwant.TranslationalSymmetry((-1, 0))
    lead_left = kwant.Builder(sym)
    lead_left[lat(0, 0)] = 1
    lead_left[lat.neighbors()] = -1
    syst.attach_lead(lead_left)
    syst.attach_lead(lead_left.reversed())

    return syst

syst = make_system(length=5).finalized()
kwant.plot(syst);

chemical_potential = 0
kwant.plotter.bands(syst.leads[0], show=False)
plt.plot([-np.pi, np.pi], [chemical_potential] * 2, 'k--')
plt.show()

#add the time dependent potential
def v(time, tau=8):
    if time < tau:
        return time / tau
    return 1
# plot the potential V(t)
times = np.linspace(0, 20)
#plt.plot(times, [v(t) for t in times])
#plt.xlabel(r'time $t$')
#plt.ylabel(r'time-dependent perturbation $V(t)$')
#plt.show()

density_operator = kwant.operator.Density(syst)

#initialize the many-body state
state = tkwant.manybody.State(syst, tmax=max(times))

densities = []
for time in times:
    state.evolve(time)
    density = state.evaluate(density_operator)
    densities.append(density)
    
densities = np.array(densities).T
for site, density in enumerate(densities):
    plt.plot(times, density, label='site {}'.format(site))
plt.xlabel(r'time $t$')
plt.ylabel(r'charge density $n$')
plt.legend()
plt.show()
