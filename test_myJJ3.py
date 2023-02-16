#  Josephson Junctions simulation with Kwant
# ===========================================================================
# Author: Kehang Zhu
# Date: 2022.03.30
# ------------------
# conductance of a Josephson-junction (Andreev reflection, superconducting gap)

from ctypes import sizeof
from sqlite3 import paramstyle
import kwant

import tinyarray
import numpy as np
# import plotly.graph_objects as goconda install matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

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
IV_B=[]

def make_system(a, W, L, mu, Delta, Metal_width, t, B):
    # a = 1;
    lat = kwant.lattice.square(norbs=2);
    syst = kwant.Builder();

    # t = 1.0;
    # W = 10;
    # L = 30;
    # mu=0.8;
    # Delta=0.2;
    # Metal_width=6;
    # B = 0.3;

    # Define the Left SC region
    for i in range(Metal_width):
        for j in range(W):
            # On-site Hamiltonian
            syst[lat(i, j)] =  (4 * t - mu) * tau_z # + Delta * tau_x

            # Hopping in y-direction
            if j > 0:
                syst[lat(i, j), lat(i, j - 1)] = - t * np.exp(1j * B * i) * tau_z

            # Hopping in x-direction
            if i > 0:
                syst[lat(i, j), lat(i - 1, j)] = -t* tau_z  
                   
    #define the normal region
    # for i in range(12, 12+Metal_width):
    #     for j in range(W):
    #         # On-site Hamiltonian
    #         syst[lat(i, j)] =  (4 * t - mu) * tau_z



    #define the left lead
    sym_left_lead = kwant.TranslationalSymmetry((-a, 0))
    left_lead = kwant.Builder(sym_left_lead,  conservation_law=-tau_z, particle_hole=tau_y)
    for j in range(W):
        left_lead[lat(0, j)] = (4 * t - mu) * tau_z #+ Delta * tau_x
        if j > 0:
            left_lead[lat(0, j), lat(0, j - 1)] = - t * tau_z
        left_lead[lat(1, j), lat(0, j)] = - t * tau_z
        
    syst.attach_lead(left_lead)

    #define the right lead
    sym_right_lead = kwant.TranslationalSymmetry((a, 0))
    right_lead = kwant.Builder(sym_right_lead, conservation_law=-tau_z, particle_hole=tau_y)

    for j in range(W):
        right_lead[lat(0, j)] = (4 * t - mu) * tau_z #+ Delta * tau_x
        if j > 0:
            right_lead[lat(0, j), lat(0, j - 1)] = -t * tau_z
        right_lead[lat(1, j), lat(0, j)] = -t * tau_z

    syst.attach_lead(right_lead)
    
    return syst


# # Finalize the system
# syst = syst.finalized()


#plot the system
# kwant.plot(syst)

# Now that we have the system, we can compute conductance
def plot_conductance(syst, energies, j):#, magnetic_field):

    data = []
    sigma = 0
    current_ =[]

    for energy in energies:
        # compute the scattering matrix at a given energy
        smatrix = kwant.smatrix(syst, energy)#  params=dict(B=Mag))
        # Conductance is N - R_ee + R_he
        sigma = smatrix.submatrix((0, 0), (0, 0)).shape[0] - smatrix.transmission((0, 0), (0, 0)) + smatrix.transmission((0, 1), (0, 0))
        # data.append(sigma)
  

        # compute the transmission probability from lead 0 to
        # lead 1
        # energies.append(energy)
        #sigma = smatrix.transmission(0, 1)
        data.append(sigma)
        current_.append(sigma * energy)
    
    IV_B.append(current_)

    # Use matplotlib to write output
    # We should see conductance steps
    # pyplot.figure()
    # pyplot.plot(energies, data)
    # pyplot.xlabel("energy [t]")
    # pyplot.ylabel("conductance [e^2/h]")
    # pyplot.show()
    # pyplot.figure()
    # pyplot.plot(energies, current_)
    # pyplot.xlabel("voltage [t]")
    # pyplot.ylabel("current")
    # pyplot.show()




def test_main():
    Mags = 0.2*np.linspace(0, 40, 40)
    energies=0.01*np.linspace(0, 50, 50)
    for j in range(0, 40):
        Bj = Mags[j];
        syst = make_system(a=1, W=10, L=30,
                mu=0.8, Delta=0.2, Metal_width=6, t=1.0, B=Bj)

        # Check that the system looks as intended.
        # kwant.plot(syst)

        # Finalize the system.
        syst = syst.finalized()

        # Check particle-hole symmetry of the scattering matrix
        # check_PHS(syst)

        # Compute and plot the conductance
        plot_conductance(syst,energies, j)
    # pyplot.figure()
    # fig = go.Figure(data=[go.Mesh3d(x=Mags, y= energies, z=IV_B, color='lightpink', opacity=0.50)])
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    Y = Mags
    X = energies
    X, Y = np.meshgrid(X, Y)
    # Plot the surface.
    # print(IV_B)
    # print(X,Y)
    # Z = np.sqrt(X**2 + Y**2)
    # print(Z)
    Z = np.array(IV_B)
    # print(Z)
    # Z = Z.flatten()
    surf = ax.plot_surface(X, Y,  Z, cmap=cm.coolwarm,  linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
    # fig = go.Figure(data= [go.Surface(x=Mags, y= energies, z=IV_B)])
     # fig.show() 
    
# Call the main function if the script gets executed (as opposed to imported).
if __name__ == '__main__':
    test_main()