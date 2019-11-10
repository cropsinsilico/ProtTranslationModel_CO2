#!/usr/bin/python
import os
import pandas
import numpy as np
from scipy import integrate

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from cis_interface.interface.CisInterface import *


_save_trajectories = False
_plot_trajectories = False


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

####Need to add lines of code to add a 3rd input file 'GRN_input.txt' and modify it's columns to create a new mRNA exp column in data_static.txt
def GrCM(grn_input, data_static, data_inp,
         save_trajectories=True, plot_trajectories=True):
    r"""Integrate input, getting trajectories for the integration.
    Args:
        data_static (pandas.DataFrame): Static information.
        data_inp (pandas.DataFrame): Initial values.
        save_trajectories (bool, optional): If True, the trajectories will
            be output to disk. Defaults to True.
        plot_trajectories (bool, optional): If True, the trajectories will
            be plotted. Defaults to True.
    Outputs:
        pandas.DataFrame: The final values after integration.
    """

    r1 = integrate.ode(Protein_translation_RNA).set_integrator(
        'vode', method='bdf')
    r2 = integrate.ode(Protein_translation_RNA).set_integrator(
        'vode', method='bdf')
    r3 = integrate.ode(Protein_translation_RNA).set_integrator(
        'vode', method='bdf')
    r4 = integrate.ode(Protein_translation_RNA).set_integrator(
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

    percentageChange = grn_input.transpose()
    percentageChange.columns = ['Glyma_ID','Averages_Elevated','Averages_ambient']

    data_static = pandas.merge(data_static, percentageChange, on='Glyma_ID')
    data_static["Averages_Elevated"] = pandas.to_numeric(data_static["Averages_Elevated"])
    data_static["Averages_ambient"] = pandas.to_numeric(data_static["Averages_ambient"])
    data_static['mRNA_eleMutant'] = data_static['mRNA_ele'] + (data_static['mRNA_ele']*data_static['Averages_Elevated']/100)
    data_static['mRNA_ambMutant'] = data_static['mRNA_Amb'] + (data_static['mRNA_Amb']*data_static['Averages_ambient']/100)
    data_static['MutantEleByWTAmb'] = data_static['Initial_protein_contentAmb']*data_static['mRNA_eleMutant']/data_static['mRNA_Amb']
    data_static['MutantAmbByWTAmb'] = data_static['Initial_protein_contentAmb']*data_static['mRNA_ambMutant']/data_static['mRNA_Amb']                         
    initial_protein_concAmb = data_static["Initial_protein_contentAmb"]
    initial_protein_concEle = data_static["Initial_protein_contentEle"]
    initial_protein_concEleMutant = data_static["MutantEleByWTAmb"]
    initial_protein_concAmbMutant = data_static["MutantAmbByWTAmb"]

    ###
    L = data_inp["L"]  # protein synthesis rate per day
    U = data_inp["U"]  # protein degradation rate per day
    # factor affecting feedback from protein concentration to rate of protein
    # synthesis from mRNA
    D = data_inp["D"]

    mRNA = data_static["Initial_protein_contentAmb"]
    r1.set_initial_value(initial_protein_concAmb,
                         t_start).set_f_params(L, U, D, mRNA)
    mRNA = data_static["Initial_protein_contentEle"]
    r2.set_initial_value(initial_protein_concEle,
                         t_start).set_f_params(L, U, D, mRNA)
    mRNA = data_static["initial_protein_concEleMutant"]
    r3.set_initial_value(initial_protein_concEleMutant,
                         t_start).set_f_params(L, U, D, mRNA)
    mRNA = data_static["initial_protein_concAmbMutant"]
    r4.set_initial_value(initial_protein_concAmbMutant,
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
        
    k = 1
    list3 = r3.y[:]
    while r3.successful() and k < num_steps:
        r3.integrate(r3.t + delta_t)
        list3 = np.vstack((list3, r3.y[:]))
        k += 1
    
    k = 1        
    list4 = r4.y[:]
    while r4.successful() and k < num_steps:
        r4.integrate(r4.t + delta_t)
        list4 = np.vstack((list4, r4.y[:]))
        k += 1

    # Make separate lists for final protein concentration from ODE results of
    # Ambient and Elevated mRNA data respectively
    #data_static.loc[:, 'new_prot_levels_Amb'] = list1[-1, :]
    #data_static.loc[:, 'new_prot_levels_Ele'] = list2[-1, :]
    data_static.loc[:, 'Ele:Amb'] = list2[-1, :] / list1[-1, :]
    data_static.loc[:, 'EleMutant:Amb'] = list3[-1, :] / list1[-1, :]
    data_static.loc[:, 'AmbMutant:Amb'] = list4[-1, :] / list1[-1, :]
    if save_trajectories:
        data_static.to_csv(
            os.path.join(TDIR, "..", "Output", "GrCM_output.txt"),
            index=False, sep="\t")

    # All done!  Plot the trajectories in separate plots:
    if plot_trajectories:
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
    return data_static

    
def main():
    # read input data as pandas df & run with output trajectories and plots
    # data_static = pandas.read_csv(
    #     os.path.join(TDIR, "..", "Input", "GrCM_static.txt"), sep="\t")
    # data_inp = pandas.read_csv(
    #     os.path.join(TDIR, "..", "Input", "GrCM_input.txt"), sep="\t")
    # GrCM(data_static, data_inp)

    # Get input from chanels (supplied by file or another model)
    in1 = CisPandasInput('GrCM_input1')
    in2 = CisPandasInput('grn_input')
    in3 = CisPandasInput('GrCM_static')
    out1 = CisPandasOutput('GrCM_output')

    flag, data_static = in3.recv()
    if not flag:
        raise Exception("GrCM: Error receiving from GrCM_static")

    while True:
        flag, data_inp = in1.recv()
        if not flag:
            print("GrCM: No more input from GrCM_input1")
        flag, grn_input = in2.recv()
        if not flag:
            print("GrCM: No more input from grn_input")
            break
        data_static = GrCM(data_static, data_inp, save_trajectories=False,
                           plot_trajectories=False)
        out1.send(data_static)

        
# The ``driver`` that will integrate the ODE(s):
if __name__ == '__main__':
    main()
