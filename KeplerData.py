#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 14:07:17 2018

@author: trustinhenderson
"""

import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def get_data():
    
    all_planets = pd.read_csv("planets_v06.csv")
    names = all_planets['name_planet']
    periods = all_planets['period_days']
    mass_types = all_planets['mass_type']
    masses = all_planets['mass_planet_jup']
    radii = all_planets['radius_planet_jup']

    restricted_period = all_planets.loc[all_planets['period_days'] <= 200]
    restricted_period_mass =  restricted_period.loc[restricted_period['mass_planet_jup'] < 0.5]

    
    return restricted_period_mass
    
"""
restricted_period_mass = get_data()


fig, obs = plt.subplots()
obs.scatter(restricted_period_mass['period_days'], restricted_period_mass['mass_planet_jup'], color='black')
obs.set(xlabel='Day', ylabel='Mass or Msini Jup')
obs.grid(True)
plt.show()
"""