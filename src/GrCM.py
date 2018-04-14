#!/usr/bin/python
import os
import sys
import csv
import numpy as np
from six import StringIO
from scipy import integrate
from matplotlib.pylab import *
import matplotlib.pyplot as plt
from pandas import *
import matplotlib.gridspec as gridspec
import signal
from cis_interface.interface.CisInterface import *
import traceback
import shutil


_save_trajectories = False
_plot_trajectories = False


def Protein_translation_Amb(t, y, data, mRNAData):
    """
    Defines ODE function Conversion of Amb_mRNA to protein
    p1,p2,p3....: Protein concentrations for all ODEs
    It will have list of all parameter values for all my ODEs, so 36 values : L, U, D for each mRNA to protein conversion equation
    
    """
    
    data = pandas.read_csv(data, sep="\t")
    
    L = data["L"].tolist() # protein synthesis rate per day
    U = data["U"].tolist() # protein degradation rate per day
    D = data["D"].tolist() # factor affecting feedback from protein
                           # concentration to rate of protein synthesis from mRNA
            
    mRNAData = pandas.read_csv(mRNAData, sep="\t")
    
    mRNA = mRNAData["mRNA_Amb"].tolist()
 
    # Output from ODE function must be a COLUMN vector, with n rows
    dydt = []
    for x in range(0,len(L)):
        temp = ((L[x]*mRNA[x])/(1+y[x]/D[x]))-(U[x]*y[x])
        dydt.append(temp)
    return dydt   

def Protein_translation_Ele(t, y, data, mRNAData):
    """
    Defines ODE function Conversion of Ele_mRNA to protein
    p1,p2,p3....: Protein concentrations for all ODEs
    It will have list of all parameter values for all my ODEs, so 36 values : L, U, D for each mRNA to protein conversion equation
    
    """
    
    data = pandas.read_csv(data, sep = "\t")
    
    L = data["L"].tolist() # protein synthesis rate per day
    U = data["U"].tolist() # protein degradation rate per day
    D = data["D"].tolist() # factor affecting feedback from protein
                           # concentration to rate of protein synthesis from mRNA
            
    mRNAData = pandas.read_csv(mRNAData, sep="\t")
    
    mRNA = mRNAData["mRNA_ele"].tolist()
 
    # Output from ODE function must be a COLUMN vector, with n rows
    dydt = []
    for x in range(0,len(L)):
        temp = ((L[x]*mRNA[x])/(1+y[x]/D[x]))-(U[x]*y[x])
        dydt.append(temp)
    return dydt  

"""
    Have the ODE functions called and have the initial values of time, end time, dt value and initial values for which we are running ODE solver.  
    # Start by specifying the integrator:
    # use ``vode`` with "backward differentiation formula"
"""

def main():
    # Setup channels
    in1 = CisPandasInput('GrCM_input1')
    in2 = CisPandasInput('GrCM_static')
    out1 = CisPandasOutput('GrCM_output')

    # Handle file names for temporary files
    dir = "./Output/Temp/"
    if os.path.exists(dir):
    	shutil.rmtree(dir)
    os.makedirs(dir)
    fname_input = os.path.join(dir, "temp_input.txt")
    fname_static = os.path.join(dir, "temp_static.txt")

    # Establish static variables
    flag, data2 = in2.recv()
    assert(flag)
    data2.to_csv(fname_static, sep="\t", index=False)
    data = data2

    # Loop over input1
    flag = True
    loop_count = 0
    while flag:

        # Get input1
        flag, data1 = in1.recv()
        if not flag:
            break
        loop_count += 1
        data1.to_csv(fname_input, sep="\t", index=False)

        # TODO: Does 'Initial_protein_content' need to be updated?

        # Set up integrator
        r1 = integrate.ode(Protein_translation_Amb).set_integrator('vode', method='bdf')
        r2 = integrate.ode(Protein_translation_Ele).set_integrator('vode', method='bdf')
        # Set the time range
        t_start = 0.0
        t_final = 200.0
        delta_t = 0.1
        # Number of time steps: 1 extra for initial condition
        num_steps = int(np.floor((t_final - t_start)/delta_t) + 1)

        # Set initial condition(s): for integrating variable and time!
        # Get the initial protein concentrations from Yu's file
        initial_protein_conc = data["Initial_protein_content"].tolist()
        r1.set_initial_value(initial_protein_conc, t_start).set_f_params(
            fname_input, fname_static)
        r2.set_initial_value(initial_protein_conc, t_start).set_f_params(
            fname_input, fname_static)

        # Additional Python step: create vectors to store trajectories
        list1 = []
        for x in range(0, len(initial_protein_conc)):
            temp = "p_" + str(x)
            vars()[temp] = np.zeros((num_steps, 1))
            list1.append(vars()[temp]) 
        list2 = []
        for x in range(0, len(initial_protein_conc)):
            temp = "p_" + str(x)
            vars()[temp] = np.zeros((num_steps, 1))
            list2.append(vars()[temp]) 
        for x in range(0, len(initial_protein_conc)):
            list1[x][0] =  initial_protein_conc[x]
        for x in range(0, len(initial_protein_conc)):
            list2[x][0] =  initial_protein_conc[x]

        # Integrate the ODE(s) across each delta_t timestep
        # First ODE
        t = np.zeros((num_steps, 1))
        t[0] = t_start
        k = 1
        while r1.successful() and k < num_steps:
            r1.integrate(r1.t + delta_t)
            t[k] = r1.t
            for x in range(0, len(initial_protein_conc)):
                list1[x][k] = r1.y[x]
            k += 1
        # Second ODE
        t = np.zeros((num_steps, 1))
        t[0] = t_start
        k = 1
        while r2.successful() and k < num_steps:
            r2.integrate(r2.t + delta_t)
            t[k] = r2.t
            for x in range(0, len(initial_protein_conc)):
                list2[x][k] = r2.y[x]
            k += 1

        # Make separate lists for final protein concentration from ODE results
        # of Ambient and Elevated mRNA data respectively
        b = []
        for j in list1:
            b.append(np.array(j[-1]).tolist())
        b = [j for i in b for j in i]
        data.loc[:,'new_prot_levels_Amb'] = b
        b1 = []
        for j in list2:
            b1.append(np.array(j[-1]).tolist())
        b1 = [j for i in b1 for j in i]
        data.loc[:,'new_prot_levels_Ele'] = b1
        data.loc[:,'Ele:Amb'] = [a/b for a,b in zip(b1,b)]
        # data.loc[:, 'Prot_ele_div_initial'] = [
        #     a / b for a, b in zip(b1, initial_protein_conc)]
        out1.send(data)

        # Clean up tree
        shutil.rmtree('./Output/Temp/')

        # Save trajectories
        if _save_trajectories:
            with open("output_amb_%d.csv" % loop_count, "w") as f:
                wr = csv.writer(f)
                wr.writerows(list1)
            with open("output_ele_%d.csv" % loop_count, "w") as f:
                wr = csv.writer(f)
                wr.writerows(list2)
            with open("Time_%d.csv" % loop_count, "wb") as f:
                wr = csv.writer(f)
                wr.writerows(t)

        # All done!  Plot the trajectories in separate plots:
        if _plot_trajectories:
            fig = figure()
            gs = gridspec.GridSpec(3,4)
            ax1 = plt.subplot2grid((3,4), (0, 0))
            ax1.plot(t, list1)
            ax1.set_xlim(t_start, t_final)
            ax1.set_xlabel('Time [minutes]')
            ax1.set_ylabel('Concentration')
            ax1.grid('on')

            ax2 = plt.subplot2grid((3,4), (0, 0))
            ax2.plot(t, list2, 'r')
            ax2.set_xlim(t_start, t_final)
            ax2.set_xlabel('Time [minutes]')
            ax2.set_ylabel('Concentration')
            ax2.grid('on')

            ax3 = plt.subplot2grid((3,4), (0, 2))
            ax3.plot(t, p_3, 'r')
            ax3.set_xlim(t_start, t_final)
            #ax3.set_xlabel('Time [minutes]')
            #ax3.set_ylabel('Concentration')
            ax3.grid('on')

            ax4 = plt.subplot2grid((3,4), (0, 3))
            ax4.plot(t, p_4, 'r')
            ax4.set_xlim(t_start, t_final)
            #ax4.set_xlabel('Time [minutes]')
            #ax4.set_ylabel('Concentration')
            ax4.grid('on')

            ax5 =  plt.subplot2grid((3,4), (2, 0))
            ax5.plot(t, p_5, 'r')
            ax5.set_xlim(t_start, t_final)
            ax5.set_xlabel('Time [minutes]')
            ax5.set_ylabel('Concentration')
            ax5.grid('on')

            ax6 = plt.subplot2grid((3,4), (2, 1))
            ax6.plot(t, p_6, 'r')
            ax6.set_xlim(t_start, t_final)
            #ax6.set_xlabel('Time [minutes]')
            #ax6.set_ylabel('Concentration')
            ax6.grid('on')

            ax7 = plt.subplot2grid((3,4), (2, 2))
            ax7.plot(t, p_7, 'r')
            ax7.set_xlim(t_start, t_final)
            #ax7.set_xlabel('Time [minutes]')
            #ax7.set_ylabel('Concentration')
            ax7.grid('on')

            ax8 = plt.subplot2grid((3,4), (2, 3))
            ax8.plot(t, p_8, 'r')
            ax8.set_xlim(t_start, t_final)
            #ax8.set_xlabel('Time [minutes]')
            #ax8.set_ylabel('Concentration')
            ax8.grid('on')


            fig.savefig('Protein_ODE_result_%d.png' % loop_count)
            
 
# The ``driver`` that will integrate the ODE(s):
if __name__ == '__main__':
    main()
