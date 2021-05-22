#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 00:32:28 2018

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

from Residuals_2 import get_residuals

### CONSTANTS
seconds_in_day = 86400
sun_mass = 1.98847542*10**30
jup_mass = 1.89818717*10**27
G = 6.67408*10**(-11)
AU_in_m = 1.495978707*10**11

Peg_b_JD, Peg_b_Residual_RV = get_residuals()


######
#---------SIMULATE SMALL PLANET-----------
######

# SEMI-MAJOR AXIS FROM PERIOD
def period_to_major_axis(period_days, mass_planet_jup, mass_star_sol):
    # INPUT PARAMETER CONVERSIONS TO ASTRO UNITS
    period_seconds = period_days * seconds_in_day  # converts days to seconds
    mass_star_kg = mass_star_sol * sun_mass
    mass_planet_kg = mass_planet_jup * jup_mass
    major_axis_meters = (((period_seconds ** 2) * G *\
                          (mass_star_kg + mass_planet_kg)) /\
                         (4 * (math.pi**2)))**(1/3)
    major_axis_AU = major_axis_meters / AU_in_m # converts to AU
    return major_axis_meters, major_axis_AU


# PERIOD FROM SEMI-MAJOR AXIS
def major_axis_to_period(major_axis_AU, mass_planet_jup, mass_star_sol):
    major_axis_meters = major_axis_AU * AU_in_m
    mass_star_kg = mass_star_sol * sun_mass
    mass_planet_kg = mass_planet_jup * jup_mass
    period_seconds = ((4*math.pi**2*major_axis_meters**3)\
    /(G*(mass_star_kg + mass_planet_kg)))**(1/2)
    period_days = period_seconds / seconds_in_day
    return period_days


# SEMI-AMPLITUDE, K
# Put math. in front of all trig functions, and math.asin(x) for arcsin
def find_rv_amp_k(period_days, mass_planet_jup, mass_star_sol):
    # INPUT PARAMETER CONVERSIONS TO ASTRO UNITS
    mass_star_kg = mass_star_sol * sun_mass
    mass_planet_kg = mass_planet_jup * jup_mass
    major_axis_meters, major_axis_AU = period_to_major_axis(period_days,
                                                            mass_planet_jup,
                                                            mass_star_sol)
    rv_amp_k = (mass_planet_kg) * ((G**(1/2)))\
                / ((mass_star_kg + mass_planet_kg) ** (1/2)\
                * (major_axis_meters) ** (1/2))
    return rv_amp_k


# MEAN ANOMALY
def find_mean_anomaly(time_variable_days, period_days):
    mean_anomaly_rad = ((2 * math.pi) * time_variable_days) / period_days
    return mean_anomaly_rad


# FINDS RV for RANGE of times
def find_rv_per_time(rv_amp_k, mean_anomaly_rad):
    # must use unitless true anomaly to avoid TypeError
    true_anomaly = mean_anomaly_rad 
    # From what I can tell, this phase shift of +0.25 is necessary entirely
    # because of astronomical convention which set phase 0 =
    # used when plotting radial velocity curves. Because I simulate a range of times, this 
    # shift is necessary in order to have the simulated RVs match with the published RVs.
    phase = true_anomaly/(2 * math.pi)
    phase = phase #+ 0.25
    true_anomaly = phase * 2 * math.pi
    rv_per_time = rv_amp_k * (math.cos(true_anomaly))
    return rv_per_time


# FINDS RV for SELECTED times and adds JITTER (51 PEG B RESIDUALS)
def find_rv_for_given_dates(dates, mass_star_sol, mass_planet_jup, period_days, e):
    output_rvs = [main(mass_star_sol, mass_planet_jup, period_days, e, date) for date in dates]
    inserted_planet_rvs = output_rvs
    output_rvs = list(map(add, output_rvs, Peg_b_Residual_RV))
    return output_rvs, inserted_planet_rvs


# def main is where all the top level running takes place
def main(mass_star_sol, mass_planet_jup, period_days, e, time_variable):

    # CONVERTS PERIOD TO SEMI-MAJOR AXIS
    major_axis_meters, major_axis_AU = period_to_major_axis(period_days,
                                                            mass_planet_jup,
                                                            mass_star_sol)

    # FINDS RV SEMI-AMPLITUDE, K
    rv_amp_k = find_rv_amp_k(period_days, mass_planet_jup, mass_star_sol)

    # FINDS MEAN ANOMALY
    mean_anomaly_rad = find_mean_anomaly(time_variable, period_days)

    # FINDS RV PER TIME
    rv_per_time = find_rv_per_time(rv_amp_k, mean_anomaly_rad)

    return rv_per_time


#major_axis_meters, major_axis_AU = period_to_major_axis(121.7285, 0.04586, 1.054)
#major_axis_meters, major_axis_AU = period_to_major_axis(4.23079, 0.48524, 1.054)
#print("Major axis: ", major_axis_AU)