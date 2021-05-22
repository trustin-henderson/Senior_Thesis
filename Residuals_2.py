#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 00:19:20 2018

@author: trustinhenderson
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

def get_residuals():

    # RV Offset Parameters for Fit
    """
    Note: These offsets are taken from the best Systemic fit. It seems that a positive offset in Systemic /subtracts/ that 
    value from all of the RVs. 
    """
    offset_ELODIE = -33251.61640
    offset_LICK6 = 14.44423877
    offset_LICK8 = 4.532548991
    offset_LICK13 = 21.70665439
    offset_HIRES = -2.378331856
    offset_HARPS = -33152.52442
    offset_APF = 0.05695235452

    ELODIE = pd.read_csv('ELODIE.csv')
    ELODIE_JD = np.array(ELODIE['BJD'])
    ELODIE_RV = np.array(ELODIE['RV'])
    ELODIE_Error = np.array(ELODIE['Error'])
    # print("Num ELODIE obs: ", len(ELODIE),
         #  "Average RV error of ELODIE data: ", np.average(ELODIE_Error))

    LICK6 = pd.read_csv('LICK6.csv')
    LICK6_JD = np.array(LICK6['BJD'])
    LICK6_RV = np.array(LICK6['RV'])
    LICK6_Error = np.array(LICK6['Error'])
    # print("Num LICK6 obs: ", len(LICK6),
          # "Average RV error of LICK6 data: ", np.average(LICK6_Error))
    
    LICK8 = pd.read_csv('LICK8.csv')
    LICK8_JD = np.array(LICK8['BJD'])
    LICK8_RV = np.array(LICK8['RV'])
    LICK8_Error = np.array(LICK8['Error'])
    # print("Num LICK8 obs: ", len(LICK8),
          # "Average RV error of LICK8 data: ", np.average(LICK8_Error))

    LICK13 = pd.read_csv('LICK13.csv')
    LICK13_JD = np.array(LICK13['BJD'])
    LICK13_RV = np.array(LICK13['RV'])
    LICK13_Error = np.array(LICK13['Error'])
    # print("Num LICK13 obs: ", len(LICK13),
          # "Average RV error of LICK13 data: ", np.average(LICK13_Error))

    HIRES = pd.read_csv('HIRES.csv')
    HIRES_JD = np.array(HIRES['BJD'])
    HIRES_RV = np.array(HIRES['RV'])
    HIRES_Error = np.array(HIRES['Error'])
    # print("Num HIRES obs: ", len(HIRES),
          # "Average RV error of HIRES data: ", np.average(HIRES_Error))

    HARPS = pd.read_csv('HARPS.csv')
    HARPS_JD = np.array(HARPS['BJD'])
    HARPS_RV = np.array(HARPS['RV'])
    HARPS_Error = np.array(HARPS['Error'])
    # print("Num HARPS obs: ", len(HARPS),
         #  "Average RV error of HARPS data: ", np.average(HARPS_Error))

    APF = pd.read_csv('51Peg_3_APF_copy_v04.csv')
    APF_JD = np.array(APF['BJD'])
    APF_RV = np.array(APF['RV'])
    APF_Error = np.array(APF['Error'])
    # print("Num APF obs: ", len(APF),
          # "Average RV error of APF data: ", np.average(APF_Error))
    
    # UPDATES RV FROM OFFSET VALUES
    # ELODIE
    new_ELODIE_RV = []
    ELODIE_decimal_precision = 0 # places after the decimal in RV obs
    for RV in ELODIE_RV:
        RV = RV - offset_ELODIE
        RV = round(RV, ELODIE_decimal_precision)  # this rounds the new, offset RVs to the nearest tenth
        RV = float(format(RV, '.0f'))  # cuts off the many extra zeros after rounding
        new_ELODIE_RV.append(RV)
    new_ELODIE_RV = np.array(new_ELODIE_RV)    

    # LICK6
    new_LICK6_RV = []
    LICK6_decimal_precision = 2 # places after the decimal in RV obs
    for RV in LICK6_RV:
        RV = RV - offset_LICK6
        RV = round(RV, LICK6_decimal_precision)  # this rounds the new, offset RVs to the nearest tenth
        RV = float(format(RV, '.2f'))  # cuts off the many extra zeros after rounding
        new_LICK6_RV.append(RV)
    new_LICK6_RV = np.array(new_LICK6_RV)    

    # LICK8
    new_LICK8_RV = []
    LICK8_decimal_precision = 2 # places after the decimal in RV obs
    for RV in LICK8_RV:
        RV = RV - offset_LICK8
        RV = round(RV, LICK8_decimal_precision)  # this rounds the new, offset RVs to the nearest tenth
        RV = float(format(RV, '.2f'))  # cuts off the many extra zeros after rounding
        new_LICK8_RV.append(RV)
    new_LICK8_RV = np.array(new_LICK8_RV)

    # LICK13
    new_LICK13_RV = []
    LICK13_decimal_precision = 2 # places after the decimal in RV obs
    for RV in LICK13_RV:
        RV = RV - offset_LICK13
        RV = round(RV, LICK13_decimal_precision)  # this rounds the new, offset RVs to the nearest tenth
        RV = float(format(RV, '.2f'))  # cuts off the many extra zeros after rounding
        new_LICK13_RV.append(RV)
    new_LICK13_RV = np.array(new_LICK13_RV)

    # HIRES
    new_HIRES_RV = []
    HIRES_decimal_precision = 8 # places after the decimal in RV obs
    for RV in HIRES_RV:
        RV = RV - offset_HIRES
        RV = round(RV, HIRES_decimal_precision)  # this rounds the new, offset RVs to the nearest tenth
        RV = float(format(RV, '.8f'))  # cuts off the many extra zeros after rounding
        new_HIRES_RV.append(RV)
    new_HIRES_RV = np.array(new_HIRES_RV)  

    # HARPS
    new_HARPS_RV = []
    HARPS_decimal_precision = 7 # places after the decimal in RV obs
    for RV in HARPS_RV:
        RV = RV - offset_HARPS
        RV = round(RV, HARPS_decimal_precision)  # this rounds the new, offset RVs to the nearest hundredth
        RV = float(format(RV, '.7f'))  # cuts off the many extra zeros after rounding, must match decimal precision
        new_HARPS_RV.append(RV)
    new_HARPS_RV = np.array(new_HARPS_RV)

    # APF
    new_APF_RV = []
    APF_decimal_precision = 2 # places after the decimal in RV obs
    for RV in APF_RV:
        RV = RV - offset_APF
        RV = round(RV, APF_decimal_precision)  # this rounds the new, offset RVs to the nearest hundredth
        RV = float(format(RV, '.2f'))  # cuts off the many extra zeros after rounding, must match decimal precision
        new_APF_RV.append(RV)
    new_APF_RV = np.array(new_APF_RV)  

    # FINDS PHASE FOR ALL OBSERVATION
    t_zero = 2440000
    period = 4.2308

    phase_ELODIE = [(((JD - t_zero)/period) -  math.floor((JD - t_zero)/period)) for JD in ELODIE_JD]

    phase_LICK6 = [(((JD - t_zero)/period) -  math.floor((JD - t_zero)/period)) for JD in LICK6_JD]

    phase_LICK8 = [(((JD - t_zero)/period) -  math.floor((JD - t_zero)/period)) for JD in LICK8_JD]

    phase_LICK13 = [(((JD - t_zero)/period) -  math.floor((JD - t_zero)/period)) for JD in LICK13_JD]

    phase_HIRES = [(((JD - t_zero)/period) -  math.floor((JD - t_zero)/period)) for JD in HIRES_JD]

    phase_HARPS = [(((JD - t_zero)/period) - math.floor((JD - t_zero)/period)) for JD in HARPS_JD]

    phase_APF = [(((JD - t_zero)/period) -  math.floor((JD - t_zero)/period)) for JD in APF_JD]

    amplitude = 55.58766012
    phaseshift = 0.645  # how do I derive this phaseshift?
    sin_x = np.linspace(-0.15, 1.15, 1000)
    sin_y = [(amplitude*np.sin(2*math.pi*x+phaseshift*period)) for x in sin_x]
    
    ELODIE_sin_y = [(amplitude*np.sin(2*math.pi*phase+phaseshift*period)) for phase in phase_ELODIE]
    residual_ELODIE_RV = []
    for i in range(len(ELODIE_sin_y)):
        residual = new_ELODIE_RV[i] - ELODIE_sin_y[i]
        residual_ELODIE_RV.append(residual)
    
    LICK6_sin_y = [(amplitude*np.sin(2*math.pi*phase+phaseshift*period)) for phase in phase_LICK6]
    residual_LICK6_RV = []
    for i in range(len(LICK6_sin_y)):
        residual = new_LICK6_RV[i] - LICK6_sin_y[i]
        residual_LICK6_RV.append(residual)
    
    LICK8_sin_y = [(amplitude*np.sin(2*math.pi*phase+phaseshift*period)) for phase in phase_LICK8]
    residual_LICK8_RV = []
    for i in range(len(LICK8_sin_y)):
        residual = new_LICK8_RV[i] - LICK8_sin_y[i]
        residual_LICK8_RV.append(residual)
    
    LICK13_sin_y = [(amplitude*np.sin(2*math.pi*phase+phaseshift*period)) for phase in phase_LICK13]
    residual_LICK13_RV = []
    for i in range(len(LICK13_sin_y)):
        residual = new_LICK13_RV[i] - LICK13_sin_y[i]
        residual_LICK13_RV.append(residual)
    
    HIRES_sin_y = [(amplitude*np.sin(2*math.pi*phase+phaseshift*period)) for phase in phase_HIRES]
    residual_HIRES_RV = []
    for i in range(len(HIRES_sin_y)):
        residual = new_HIRES_RV[i] - HIRES_sin_y[i]
        residual_HIRES_RV.append(residual)
    
    APF_sin_y = [(amplitude*np.sin(2*math.pi*phase+phaseshift*period)) for phase in phase_APF]
    residual_APF_RV = []
    for i in range(len(APF_sin_y)):
        residual = new_APF_RV[i] - APF_sin_y[i]
        residual_APF_RV.append(residual)
    
    HARPS_sin_y = [(amplitude*np.sin(2*math.pi*phase+phaseshift*period)) for phase in phase_HARPS]
    residual_HARPS_RV = []
    for i in range(len(HARPS_sin_y)):
        residual = new_HARPS_RV[i] - HARPS_sin_y[i]
        residual_HARPS_RV.append(residual)
    
    
    all_residual_JD = np.concatenate((ELODIE_JD, LICK6_JD, LICK8_JD, LICK13_JD, HIRES_JD, HARPS_JD, APF_JD), axis=0)
    # print(Peg_b_JD)
    all_residual_RV = np.concatenate((residual_ELODIE_RV, residual_LICK6_RV, residual_LICK8_RV, 
                                      residual_LICK13_RV, residual_HIRES_RV, residual_HARPS_RV,
                                      residual_HARPS_RV, residual_APF_RV), axis=0)
    
    all_JD = all_residual_JD
    all_RV = all_residual_RV  # data
    
    alltimeorder = np.argsort(all_JD)
    # print(timeorder)
    all_residual_JD = all_JD[alltimeorder]
    all_residual_RV = all_RV[alltimeorder]


    return all_residual_JD, all_residual_RV