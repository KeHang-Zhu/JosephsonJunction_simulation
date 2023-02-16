# Tutorial 2.4.2. Closed systems
# ==============================
#
# Physics background
# ------------------
#  Fock-darwin spectrum of a quantum dot (energy spectrum in
#  as a function of a magnetic field)
#
# Kwant features highlighted
# --------------------------
#  - Use of `hamiltonian_submatrix` in order to obtain a Hamiltonian
#    matrix.

from cmath import exp
import numpy as np
import kwant

# For eigenvalue computation
import scipy.sparse.linalg as sla

# For plotting
from matplotlib import pyplot


def make_system(a=1, t=1.0, r=10):
    # Start with an empty tight-binding system and a single square lattice.
    # `a` is the lattice constant (by default set to 1 for simplicity).

    lat = kwant.lattice.square(a, norbs=1)

    syst = kwant.Builder()

    # Define the quantum dot
    def circle(pos):
        (x, y) = pos
        rsq = x ** 2 + y ** 2
        return rsq < r ** 2

    def hopx(site1, site2, B):
        # The magnetic field is controlled by the parameter B
        y = site1.pos[1]
        return -t * exp(-1j * B * y)

    syst[lat.shape(circle, (0, 0))] = 4 * t
    # hoppings in x-direction
    syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopx
    # hoppings in y-directions
    syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = -t

    # It's a closed system for a change, so no leads
    return syst


def plot_spectrum(syst, Bfields):

    energies = []
    for B in Bfields:
        # Obtain the Hamiltonian as a sparse matrix
        ham_mat = syst.hamiltonian_submatrix(params=dict(B=B), sparse=True)

        # we only calculate the 15 lowest eigenvalues
        ev = sla.eigsh(ham_mat.tocsc(), k=15, sigma=0,
                       return_eigenvectors=False)

        energies.append(ev)

    pyplot.figure()
    pyplot.plot(Bfields, energies)
    pyplot.xlabel("magnetic field [arbitrary units]")
    pyplot.ylabel("energy [t]")
    pyplot.show()

def sorted_eigs(ev):
    evals, evecs = ev
    evals, evecs = map(np.array, zip(*sorted(zip(evals, evecs.transpose()))))
    return evals, evecs.transpose()

def plot_wave_function(syst, B=0.001):
    # Calculate the wave functions in the system.
    ham_mat = syst.hamiltonian_submatrix(sparse=True, params=dict(B=B))
    evals, evecs = sorted_eigs(sla.eigsh(ham_mat.tocsc(), k=20, sigma=0))

    # Plot the probability density of the 10th eigenmode.
    kwant.plotter.map(syst, np.abs(evecs[:, 9])**2,
                      colorbar=False, oversampling=1)


def plot_current(syst, B=0.001):
    # Calculate the wave functions in the system.
    ham_mat = syst.hamiltonian_submatrix(sparse=True, params=dict(B=B))
    evals, evecs = sorted_eigs(sla.eigsh(ham_mat.tocsc(), k=20, sigma=0))

    # Calculate and plot the local current of the 10th eigenmode.
    J = kwant.operator.Current(syst)
    current = J(evecs[:, 9], params=dict(B=B))
    kwant.plotter.current(syst, current, colorbar=False)


def test_main():
    syst = make_system()

    # Check that the system looks as intended.
    #kwant.plot(syst)

    # Finalize the system.
    syst = syst.finalized()

    # The following try-clause can be removed once SciPy 0.9 becomes uncommon.
    try:
        # We should observe energy levels that flow towards Landay
        # level energies with increasing magnetic field.
        plot_spectrum(syst, [iB * 0.002 for iB in range(100)])

        # Plot an eigenmode of a circular dot. Here we create a larger system for
        # better spatial resolution.
        syst = make_system(r=30).finalized()
        plot_wave_function(syst)
        plot_current(syst)
    except ValueError as e:
        if e.message == "Input matrix is not real-valued.":
            print("The calculation of eigenvalues failed because of a bug in SciPy 0.9.")
            print("Please upgrade to a newer version of SciPy.")
        else:
            raise


# Call the main function if the script gets executed (as opposed to imported).
# See <http://docs.python.org/library/__main__.html>.
if __name__ == '__main__':
    test_main()
