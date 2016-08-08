#!/usr/bin/python
import os
import pandas
import numpy as np
from scipy import integrate

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def Protein_translation_RNA(t, y, L, U, D, mRNA):
    """
    Defines ODE function Conversion of Amb_mRNA to protein
    p1,p2,p3....: Protein concentrations for all ODEs
    It will have list of all parameter values for all my ODEs, so 36 values :
    L, U, D for each mRNA to protein conversion equation

    """
    # Output from ODE function must be a COLUMN vector, with n rows
    return (L * mRNA) / (1.0 + y / D) - U * y


"""
    Have the ODE functions called and have the initial values of time, end
    time, dt value and initial values for which we are running ODE solver.
    # Start by specifying the integrator: use ``vode`` with "backward
    # differentiation formula"
"""

TDIR = os.path.dirname(os.path.abspath(__file__))


def main():
    # read input data as pandas df
    data_static = pandas.read_csv(
        os.path.join(TDIR, "..", "Input", "GrCM_static.txt"), sep="\t")
    data_inp = pandas.read_csv(
        os.path.join(TDIR, "..", "Input", "GrCM_input.txt"), sep="\t")

    r1 = integrate.ode(Protein_translation_RNA).set_integrator(
        'vode', method='bdf')
    r2 = integrate.ode(Protein_translation_RNA).set_integrator(
        'vode', method='bdf')
    # Set the time range
    t_start = 0.0
    t_final = 200.0
    delta_t = 0.1
    t = np.arange(t_start, t_final + delta_t, delta_t)
    # Number of time steps: 1 extra for initial condition
    num_steps = t.shape[0]

    # Set initial condition(s): for integrating variable and time! Get the
    # initial protein concentrations from Yu's file
    initial_protein_conc = data_static["Initial_protein_content"]

    ###
    L = data_inp["L"]  # protein synthesis rate per day
    U = data_inp["U"]  # protein degradation rate per day
    # factor affecting feedback from protein concentration to rate of protein
    # synthesis from mRNA
    D = data_inp["D"]

    mRNA = data_static["mRNA_Amb"]
    r1.set_initial_value(initial_protein_conc,
                         t_start).set_f_params(L, U, D, mRNA)
    mRNA = data_static["mRNA_ele"]
    r2.set_initial_value(initial_protein_conc,
                         t_start).set_f_params(L, U, D, mRNA)

    # Integrate the ODE(s) across each delta_t timestep
    k = 1
    list1 = r1.y[:]
    while r1.successful() and k < num_steps:
        r1.integrate(r1.t + delta_t)
        list1 = np.vstack((list1, r1.y[:]))
        k += 1

    k = 1
    list2 = r2.y[:]
    while r2.successful() and k < num_steps:
        r2.integrate(r2.t + delta_t)
        list2 = np.vstack((list2, r2.y[:]))
        k += 1

    # Make separate lists for final protein concentration from ODE results of
    # Ambient and Elevated mRNA data respectively
    data_static.loc[:, 'new_prot_levels_Amb'] = list1[-1, :]
    data_static.loc[:, 'new_prot_levels_Ele'] = list2[-1, :]
    data_static.loc[:, 'Ele:Amb'] = list2[-1, :] / list1[-1, :]
    data_static.to_csv(
        os.path.join(TDIR, "..", "Output", "GrCM_output.txt"),
        index=False, sep="\t")

    # All done!  Plot the trajectories in separate plots:
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 4)
    ax1 = fig.add_subplot(gs[0])
    ax1.semilogx(t, list1)
    ax1.set_xlim(t_start, t_final)
    ax1.set_xlabel('Time [minutes]')
    ax1.set_ylabel('Concentration')
    ax1.grid('on')

    ax2 = fig.add_subplot(gs[1])
    ax2.semilogx(t, list2, 'r')
    ax2.set_xlim(t_start, t_final)
    ax2.set_xlabel('Time [minutes]')
    ax2.set_ylabel('Concentration')
    ax2.grid('on')

    for i in range(2, 8):
        ax = fig.add_subplot(gs[i])
        ax.semilogx(t, list2[:, i + 1], 'r')
        ax.set_xlim(t_start, t_final)
        ax.set_xlabel('Time [minutes]')
        ax.grid('on')

    plt.savefig('Protein_ODE_result.png', bbox='tight')

# The ``driver`` that will integrate the ODE(s):
if __name__ == '__main__':
    main()
