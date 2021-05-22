#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 02:36:21 2018

@author: trustinhenderson
"""

import time

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import constants as const
from astropy import units as u
from astropy.stats import LombScargle as ls
import math
from operator import add
import matplotlib as mpl
from matplotlib import cm as CM
from matplotlib import mlab as ML
import multiprocessing as mp

from Residuals_2 import get_residuals
import RV_Planet_Creation_3 as RVFunc
import KeplerData


def get_periodogram(mass_planet_jup, period_days):
    # INPUT PARAMETERS FOR SYSTEM
    # mass_star_sol = float(input("Enter the star mass in solar masses: "))
    mass_star_sol = 1.054  # mass of 51 Peg in solar masses
    e = 0  # eccentricity
    dates = Peg_b_JD
        
    output_rvs, inserted_planet_rvs = RVFunc.find_rv_for_given_dates(dates,
                                                                     mass_star_sol,
                                                                     mass_planet_jup,
                                                                     period_days,
                                                                     e)
    
    # Creates and draws output plot, and gives amplitude of curve
    # print("Amplitude of RV Curve: ", rv_amp_k)
    """
    fig1, obs = plt.subplots()
    obs.scatter(dates[0:20], output_rvs[0:20], color='black')
    obs.scatter(dates[0:20], inserted_planet_rvs[0:20], color='red')
    obs.set(xlabel='Day', ylabel='RV [m/s]', title='RV per Time')
    obs.grid(True)
    obs.legend(loc='upper right')
    plt.show()
    """
    
    x = np.array(dates)
    y = np.array(output_rvs)  # data
    # print("Dates: ", type(x), x)
    # print("RVs: ", type(y), y)
    timeorder = np.argsort(x)
    # print(timeorder)
    x = x[timeorder]
    y = y[timeorder]
    dates = x
    output_rvs = y
    # print("Dates: ", dates)
    # print("RVs: ", output_rvs)
    
    t = dates
    y = output_rvs
    frequency, power = ls(t, y).autopower(minimum_frequency= 0.004, 
                         maximum_frequency=0.5,
                         samples_per_peak = 10,
                         method='fast')
    
    # the minimum frequency is the max period
    # the max frequency is the minimum period
    
    """
    fig3, lomb = plt.subplots()
    lomb.plot((1/frequency), power, color='black')
    lomb.set(xlabel='Period', ylabel='Lomb-Scargle Power', title='Periodogram')
    lomb.grid(True)
    lomb.legend(loc='upper right')
    plt.show()
    """
    
    list_frequency = [1/period_days]
    power = (ls(t, y).power(list_frequency))
    # FAP = (ls(t,y).false_alarm_probability(power))
    
    return power


def mass_period_creation(N, gridsize):
    # jupiter masses, first input MUST BE POSITIVE
    mass_range_arr = np.array(np.linspace(0.01, 0.5, N))
    
    lin_mass_arr = [q for q in mass_range_arr]
    lin_mass_arr = np.repeat(lin_mass_arr, N)
    new_mass_arr = [s for s in lin_mass_arr]

    period_range_arr = np.linspace(2, 200, N)  # days
    new_period_arr = [r for r in period_range_arr] * N

    return new_mass_arr, new_period_arr


def get_power(x):
    mass_planet_jup = new_mass_arr[x]
    period_days = new_period_arr[x]
    power = get_periodogram(mass_planet_jup, period_days)
    
    return power


if __name__ == '__main__':
    start = time.time()
    
    Peg_b_JD, Peg_b_Residual_RV = get_residuals()
    
    N = 1500
    gridsize = 150
    
    new_mass_arr, new_period_arr = mass_period_creation(N, gridsize)
    
    all_FAP_arr = []  # array where each element is the FAP at period and mass
    all_Power_arr = []
    all_points_arr = []
    
    cores = 4
    pool = mp.Pool(processes=cores)  # number of cores running
    all_points_arr = pool.map(get_power, range(0, len(new_mass_arr)))
    # print(results)

    # if 'bins=None', then color of each hexagon corresponds directly to its count
    # 'C' is optional--it maps values to x-y coordinates; if 'C' is None (default)
    # the result is a pure 2D histogram
    new_period_arr = np.asarray(new_period_arr)
    new_mass_arr = np.asarray(new_mass_arr)

    x = new_period_arr
    y = new_mass_arr
    
    restricted_period_mass = KeplerData.get_data()
    #print(restricted_period_mass['period_days'])
    #print(restricted_period_mass['mass_planet_jup'])
    
    fig5, hexplot = plt.subplots()
    graph = hexplot.hexbin(x, y, C=all_points_arr, gridsize=gridsize,
                           cmap='cool', bins=None, yscale='log') #, xscale='log', yscale='log')
    hexplot.set_yscale('log')
    hexplot.scatter(restricted_period_mass['period_days'],
                    restricted_period_mass['mass_planet_jup'],
                    color='black', alpha=0.5, s=2)
    hexplot.set_xlabel('Period in Days')
    hexplot.set_ylabel('Msini in Jupiter Masses')
    # hexplot.set_yscale('log')
    hexplot.axis([x.min(), x.max(), y.min(), y.max()]) 
    cb = fig5.colorbar(graph) #.set_label('Power')
    cb.set_label('Power')

    print("Gridsize: ", gridsize)
    print("N: ", N)
    print("# of cores: ", cores)
    
    end = time.time()
    time_to_run = (end - start) 
    print(time_to_run / 60, "minutes")
    
    points_per_sec = N**2 / time_to_run
    print(round(points_per_sec, 2), " points / second")

    fig5.savefig('plot_1500_150_4.pdf')
    
    text_file = open("plot_1500_150_4.txt", "w")
    line1 = str(time_to_run / 60) + " minutes"
    line2 = str(round(points_per_sec, 2)) + " points / second"
    line3 = "N: " + str(N)
    line4 = "Gridsize: " + str(gridsize)
    line5 = "# of cores: " + str(cores)
    text_file.write("%s \n %s \n %s \n %s \n %s \n" % (line1, line2, line3, line4, line5))
    text_file.close()
    