#!/usr/bin/env python3

import kwant
from math import sin

lat = kwant.lattice.square(a=1, norbs=1)
syst = kwant.Builder()

syst[(lat(x, 0) for x in range(20))] = 1
syst[lat.neighbors()] = -1

def onsite(site, fac, omega, time):
    return 1 + fac * sin(omega * time)

# add the time-dependent onsite potential
syst[lat(10, 0)] = onsite

lead = kwant.Builder(kwant.TranslationalSymmetry((-1, 0)))
lead[lat(0, 0)] = 1             #dont
lead[lat.neighbors()] = -1   #dont understand
syst.attach_lead(lead)
syst.attach_lead(lead.reversed())

# plot the system
kwant.plot(syst, site_color=lambda s: 'r' if s in [lat(10, 0)] else 'k',
           lead_color='grey');

syst = syst.finalized()