# Tutorial 2.3.3. Nontrivial shapes
# =================================
#
# Physics background
# ------------------
#  Flux-dependent transmission through a quantum ring
#
# Kwant features highlighted
# --------------------------
#  - More complex shapes with lattices
#  - Allows for discussion of subtleties of `attach_lead` (not in the
#    example, but in the tutorial main text)
#  - Modifcations of hoppings/sites after they have been added

from cmath import exp
from math import pi

import kwant

# For plotting
from matplotlib import pyplot


def make_system(a=1, t=1.0, W=10, r1=10, r2=20):
    # Start with an empty tight-binding system and a single square lattice.
    # `a` is the lattice constant (by default set to 1 for simplicity).

    lat = kwant.lattice.square(a)

    syst = kwant.Builder()

    #### Define the scattering region. ####
    # Now, we aim for a more complex shape, namely a ring (or annulus)
    def ring(pos):
        (x, y) = pos
        rsq = x ** 2 + y ** 2
        return (r1 ** 2 < rsq < r2 ** 2)

    # and add the corresponding lattice points using the `shape`-function
    syst[lat.shape(ring, (0, r1 + 1))] = 4 * t
    syst[lat.neighbors()] = -t

    # In order to introduce a flux through the ring, we introduce a phase on
    # the hoppings on the line cut through one of the arms.  Since we want to
    # change the flux without modifying the Builder instance repeatedly, we
    # define the modified hoppings as a function that takes the flux as its
    # parameter phi.
    def hopping_phase(site1, site2, phi):
        return -t * exp(1j * phi)

    def crosses_branchcut(hop):
        ix0, iy0 = hop[0].tag

        # builder.HoppingKind with the argument (1, 0) below
        # returns hoppings ordered as ((i+1, j), (i, j))
        return iy0 < 0 and ix0 == 1  # ix1 == 0 then implied

    # Modify only those hopings in x-direction that cross the branch cut
    def hops_across_cut(syst):
        for hop in kwant.builder.HoppingKind((1, 0), lat, lat)(syst):
            if crosses_branchcut(hop):
                yield hop
    syst[hops_across_cut] = hopping_phase

    #### Define the leads. ####
    # left lead
    sym_lead = kwant.TranslationalSymmetry((-a, 0))
    lead = kwant.Builder(sym_lead)

    def lead_shape(pos):
        (x, y) = pos
        return (-W / 2 < y < W / 2)

    lead[lat.shape(lead_shape, (0, 0))] = 4 * t
    lead[lat.neighbors()] = -t

    #### Attach the leads and return the system. ####
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst


def plot_conductance(syst, energy, fluxes):
    # compute conductance

    normalized_fluxes = [flux / (2 * pi) for flux in fluxes]
    data = []
    for flux in fluxes:
        smatrix = kwant.smatrix(syst, energy, params=dict(phi=flux))
        data.append(smatrix.transmission(1, 0))

    pyplot.figure()
    pyplot.plot(normalized_fluxes, data)
    pyplot.xlabel("flux [flux quantum]")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()


def test_main():
    syst = make_system()

    # Check that the system looks as intended.
    kwant.plot(syst)

    # Finalize the system.
    syst = syst.finalized()

    # We should see a conductance that is periodic with the flux quantum
    plot_conductance(syst, energy=0.15, fluxes=[0.01 * i * 3 * 2 * pi
                                                for i in range(100)])


# Call the main function if the script gets executed (as opposed to imported).
# See <http://docs.python.org/library/__main__.html>.
if __name__ == '__main__':
    test_main()
