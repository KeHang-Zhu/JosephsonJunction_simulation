#  Josephson Junctions simulation with Kwant
# ===========================================================================
# Author: Kehang Zhu
# Date: 2022.03.23
# ------------------
# conductance of a Josephson-junction (Andreev reflection, superconducting gap)

from sqlite3 import paramstyle
import kwant

import tinyarray
import numpy as np

# For plotting
from matplotlib import pyplot

tau_x = tinyarray.array([[0, 1], [1, 0]])
tau_y = tinyarray.array([[0, -1j], [1j, 0]])
tau_z = tinyarray.array([[1, 0], [0, -1]])

################
#  SystemInfo  #
################
#a: lattice size
#L: the width of the metal inside the JJ
#t : hopping rate inside the SC
#t2: hopping rate inside the normal metal

def make_system(a=1, W=10, L=10, barrier=1.5, t2 = 1.0,
                mu=0.4, Delta=0.1, t=1.0, B = 0.1):
    # Start with an empty tight-binding system. On each site, there
    # are now electron and hole orbitals, so we must specify the
    # number of orbitals per site. The orbital structure is the same
    # as in the Hamiltonian.
    lat = kwant.lattice.square(norbs=2)
    syst = kwant.Builder()

    #### Define the Normal metal region. ####
    syst[(lat(x, y) for x in range(L)
         for y in range(W))] = (4 * t - mu) * tau_z

    # Hoppings
    syst[lat.neighbors()] = -t * tau_z
    
    #### implement Pierce construction ####
    ##In an external uniform magnetic field applied in the z-direction, the hopping integral acquires the Peierls phase factor as ##
    
    # Modify only those hopings in x-direction that cross the branch cut
    def hopping_phase(site):
        (x1, y1) = site.family_a.pos;
        (x2, y2) = site.family_b.pos;
        if (x1 == x2):
            return -t * exp(1j * B * x1) * tau_z
        else:
            return -t * tau_z
    def hops_across_cut(syst):
        for hop in kwant.builder.HoppingKind((1, 0), lat, lat)(syst):
            syst[hop] = hopping_phase(hop)

    
    #### Define the leads. ####
    # Left lead - SC, so the order parameter is non-zero.
    sym_left = kwant.TranslationalSymmetry((-a, 0))
    # Specify the conservation law used to treat electrons and holes separately.
    # We only do this in the left lead, where the pairing is non-zero.
    lead0 = kwant.Builder(sym_left, conservation_law=-tau_z, particle_hole=tau_y)
    lead0[(lat(0, j) for j in range(W))] = (4 * t - mu) * tau_z + Delta * tau_x
    lead0[lat.neighbors()] = -t * tau_z
    # Right lead - SC
    sym_right = kwant.TranslationalSymmetry((a, 0))
    lead1 = kwant.Builder(sym_right)
    lead1[(lat(0, j) for j in range(W))] = (4 * t - mu) * tau_z + Delta * tau_x
    lead1[lat.neighbors()] = -t * tau_z

    #### Attach the leads and return the system. ####
    syst.attach_lead(lead0)
    syst.attach_lead(lead1)

    return syst


def plot_conductance(syst, energies):#, magnetic_field):
    # Compute conductance
    #normalized_B = [Mag / (2 * pi) for Mag in Mags]
    data = []
    sigma = 0;
    current_ =[]
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)#  params=dict(B=Mag))
        # Conductance is N - R_ee + R_he
        sigma = smatrix.submatrix((0, 0), (0, 0)).shape[0] - smatrix.transmission((0, 0), (0, 0)) + smatrix.transmission((0, 1), (0, 0));
        data.append(sigma)
        current_.append(sigma * energy)
    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("energy [t]")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()
    pyplot.figure()
    pyplot.plot(energies, current_)
    pyplot.xlabel("voltage [t]")
    pyplot.ylabel("current")
    pyplot.show()
    
def check_PHS(syst):
    # Scattering matrix
    s = kwant.smatrix(syst, energy=0)
    # Electron to electron block
    s_ee = s.submatrix((0,0), (0,0))
    # Hole to hole block
    s_hh = s.submatrix((0,1), (0,1))
    print('s_ee: \n', np.round(s_ee, 3))
    print('s_hh: \n', np.round(s_hh[::-1, ::-1], 3))
    print('s_ee - s_hh^*: \n',
          np.round(s_ee - s_hh[::-1, ::-1].conj(), 3), '\n')
    # Electron to hole block
    s_he = s.submatrix((0,1), (0,0))
    # Hole to electron block
    s_eh = s.submatrix((0,0), (0,1))
    print('s_he: \n', np.round(s_he, 3))
    print('s_eh: \n', np.round(s_eh[::-1, ::-1], 3))
    print('s_he + s_eh^*: \n',
          np.round(s_he + s_eh[::-1, ::-1].conj(), 3))

#######################
# Current Calculation #
#######################
# def plot_currents(syst, currents):
#     fig, axes = plt.subplots(1, len(currents))
#     if not hasattr(axes, '__len__'):
#         axes = (axes,)
#     for ax, (title, current) in zip(axes, currents):
#         kwant.plotter.current(syst, current, ax=ax, colorbar=False)
#         ax.set_title(title)
#     plt.show()


def test_main():
    syst = make_system(W=10)

    # Check that the system looks as intended.
    kwant.plot(syst)

    # Finalize the system.
    syst = syst.finalized()

    # Check particle-hole symmetry of the scattering matrix
    check_PHS(syst)

    # Compute and plot the conductance
    plot_conductance(syst, energies=[0.002 * i for i in range(-100, 100)])


# Call the main function if the script gets executed (as opposed to imported).
if __name__ == '__main__':
    test_main()
