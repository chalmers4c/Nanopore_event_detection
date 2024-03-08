# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 00:18:17 2023

The code will analyse an .abf file and detect translocation event based on simple threshold setting.
Input: .abf file
Output: 
    Translocation Events' detection and data summary
    Events overlay

@author: Chalmers Chau @ University of Leeds
"""
#%% Packages to import
"""Essential packages"""
import tkinter as tk
from tkinter import filedialog
import pyabf
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from pybaselines import Baseline, utils
import seaborn as sns
import math
import matplotlib.ticker as ticker
import scipy
from scipy import signal
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import random
import os
import glob
import pathlib

"""Optional packages"""
from palettable.wesanderson import GrandBudapest5_5
from palettable.cartocolors.qualitative import Vivid_10
from palettable.wesanderson import Mendl_4
#%% Global variables
# Trace region selection, if you have a particular region in seconds that you want to process, turn ROI = 'Yes', otherwise use 'No'

#ROI = 'Yes'
ROI_range_from = 4, #Control the region in seconds, irrelevant but is activated, because lazy.
ROI_range_to = 9

ROI = 'No'

# Standard deviation setting of the threshold.
Threshold_std = 5 # this depends on data set
Threshold_prominence = 5 # this has the same value as std
Threshold_width = 0.001 # width of the feature to be considered

# Below to control the width of the translocation event overlay, adjust if more time before and after is needed to visualise the event.
#Before_peak_time, default is 10
Before_peak_time_value = 10

#After_peak_time, default is 30
After_peak_time_value = 70

#%% Functions to load file and export files
def file_path(): 
    root = tk.Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    file_full_path = filedialog.askopenfilename()
    root.destroy()
    folder_path = os.path.dirname(file_full_path)
    return file_full_path, folder_path
file_path, folder_path = file_path()

def save_to_path_png(Path, Name): #Save plot to same path as samples
    return plt.savefig(os.path.join(Path, Name)+'.png', format = 'png', dpi = 300)

def save_to_path_svg(Path, Name):
    return plt.savefig(os.path.join(Path, Name)+'.svg', format = 'svg')

def excel_export(Path, Name):
    excel_path = os.path.join(folder_path, Name)+'.xlsx'
    with pd.ExcelWriter(excel_path) as excel_writer:
        df_Peak_search_result.to_excel(excel_writer, index = False, sheet_name = 'Events Summary')
    return
File_name = str(os.path.basename(file_path).split('.')[0])

#%% Read abf trace with pyabf
def pyabf_trace(Path):
    trace_abf = pyabf.ABF(Path)
    trace_current = trace_abf.setSweep(sweepNumber = 0, channel = 0)
    trace_current_time = trace_abf.sweepX
    trace_current_current = trace_abf.sweepY
    return trace_current_time, trace_current_current, trace_abf
Time_axis, Current_axis, trace_abf = pyabf_trace(Path = file_path)


#%% baseline tracing through baseline_fitter using Statistics-sensitive Non-linear Iterative Peak-clipping (SNIP).
def baseline_adjustment(ROI_range_from, 
                        ROI_range_to,
                        ROI,
                        max_window = 50,
                        smooth_window = 500,
                        order = 4,
                        ):    
    baseline_fitter = Baseline(Time_axis, check_finite=False)
    baseline = baseline_fitter.snip(data = Current_axis, 
                                    max_half_window = max_window, 
                                    decreasing = True, 
                                    smooth_half_window = smooth_window,
                                    filter_order = order
                                    )[0]
    baseline_corrected = Current_axis - baseline
    df = pd.DataFrame({'Time_axis': Time_axis, 
                       'Current_baseline_adjusted': baseline_corrected})
    
    if ROI == 'Yes':
        df_rows_selected = df[(df['Time_axis'] >= ROI_range_from) & (df['Time_axis'] <= ROI_range_to)]
        df_rows_selected.reset_index(drop = True, inplace = True)
        Time_min = df_rows_selected.iloc[:, 0].min()
        df_rows_selected['Time_axis'] = df_rows_selected['Time_axis'] - Time_min
    elif ROI == 'No':
        df_rows_selected = df
    else:
        print('Either Yes or No for the trace ROI')
    return baseline, baseline_corrected, df, df_rows_selected

baseline, baseline_corrected, df_baseline_corrected, df_baseline_corrected_ROI = baseline_adjustment(ROI = ROI, ROI_range_from = ROI_range_from, ROI_range_to = ROI_range_to, max_window = 50, smooth_window = 500, order = 4)
#%% Peaks_search with scipy
def Peaks_search(ROI,
                 Before_peak_time = Before_peak_time_value,
                 After_peak_time = After_peak_time_value, #default value
                 Threshold_std = Threshold_std,
                 Threshold_prominence = Threshold_prominence,
                 Threshold_width = Threshold_width):
    if ROI == 'No':
        peaks_signal = df_baseline_corrected.iloc[:, 1]  # trace after baseline adjustment
    elif ROI == 'Yes':
        peaks_signal = df_baseline_corrected_ROI.iloc[:, 1]
    else:
        print('Either Yes or No for the trace ROI')
    
    Sample_frequency = trace_abf.sampleRate # Sampling frequency in Hertz
    
    # Event window parameters the value input here will sum to the total window frame
    pretrigger_window = (Before_peak_time * Sample_frequency)/10000  # Pre-event time window in us
    posttrigger_window = (After_peak_time * Sample_frequency)/10000  # Post-event time window in us
    
    # Set parameters of the Find peaks function
    thresh_min = Threshold_std * np.std(peaks_signal) #7 standard deivation from the mean
    thresh_max = float('inf')
    thresh_prominence = Threshold_prominence * np.std(peaks_signal)
    thresh_min_width = Threshold_width * (Sample_frequency/10000)
     
    peaks, peaks_dict = find_peaks(x = -peaks_signal, 
                                   height = (thresh_min, thresh_max),  # Min and max thresholds to detect peaks.
                                   threshold = None,  # Min and max vertical distance to neighboring samples.
                                   distance = None,  # Min horizontal distance between peaks.
                                   prominence = thresh_prominence,  # Vertical distance between the peak and lowest contour line.
                                   width = thresh_min_width,  # Min required width (in bins). E.g. For 10Khz, 10 bins = 1 ms.
                                   wlen = None,  # Window length to calculate prominence.
                                   rel_height = 0.5,  # Relative height at which the peak width is measured.
                                   plateau_size = None)
         
    # Create table with results
    df_results = pd.DataFrame(columns = ['event', 
                                         'peak_position', 
                                         'peak_position_s',
                                         'event_start', 
                                         'event_end',
                                         'Peak_Amp_pA',
                                         'Width_ms',
                                         'Peak_Amp_nA',
                                         'Area_nA·ms',
                                         'Area_pA·ms'])
      
    df_results.event = np.arange(1, len(peaks) + 1)
    df_results.peak_position = peaks
    df_results.peak_position_s = peaks / Sample_frequency  # Divided by fs to get s
    df_results.event_start = peaks_dict['left_ips'] - pretrigger_window
    df_results.event_end = peaks_dict['right_ips'] + posttrigger_window
    df_results.Width_ms = peaks_dict['widths']/(Sample_frequency/1000) # Width (ms) at half-height
    df_results.Peak_Amp_pA = peaks_dict['peak_heights']  # height parameter is needed
    df_results.Peak_Amp_nA = df_results.Peak_Amp_pA/1000  # height parameter is needed
    for i, event in df_results.iterrows():
        individual_event = peaks_signal.iloc[int(event.event_start) : int(event.event_end)]
        df_results.loc[i, 'Area_pA·ms'] = -np.round(individual_event.sum(), 1)/(Sample_frequency/1000) 
    for i, event in df_results.iterrows():
        individual_event = peaks_signal.iloc[int(event.event_start) : int(event.event_end)]
        df_results.loc[i, 'Area_nA·ms'] = -np.round(individual_event.sum(), 1)/(Sample_frequency/1000)/1000
    return df_results

df_Peak_search_result = Peaks_search(ROI = ROI)

#%% Translocation peaks position data
def Peaks_data(ROI):
    # Construct dfs, one for start end time after peak tracking, one for the actual scatter of the trace
    if ROI == 'Yes':
        df_trace = pd.DataFrame({'Time' : df_baseline_corrected_ROI.iloc[:, 0], 'Current': df_baseline_corrected_ROI.iloc[:, 1]})
    elif ROI == 'No':
        df_trace = pd.DataFrame({'Time' : df_baseline_corrected.iloc[:, 0], 'Current': df_baseline_corrected.iloc[:, 1]})
    else:
        print('Either Yes or No for the trace ROI')
    
    df_start_end = pd.DataFrame({'Start': df_Peak_search_result['event_start']/100000, 
                                 'End': df_Peak_search_result['event_end']/100000})
    
    # Create empty list to store DataFrames
    dfs_peaks = {}

    # Iterate over each row in the DataFrame
    for index, row in df_start_end.iterrows():
        Event_Start = row['Start']
        Event_End = row['End']
        
        # Select data points between start and end on the x-axis
        mask = (df_trace['Time'] >= Event_Start) & (df_trace['Time'] <= Event_End)
        
        # Create a DataFrame for the current region
        df_peak = df_trace.loc[mask].copy()
        df_peak.reset_index(drop=True, inplace=True)
        dfs_peaks[f'peak{index+1}'] = df_peak
        
        # Iterate over the peak DataFrames
    for key, df in dfs_peaks.items():
        # Get the smallest value in the "Time" column
        min_time = df['Time'].min()
        # Recalculate the "Time" column relative to the smallest value
        df['Time'] = df['Time'] - min_time
    return dfs_peaks

dfs_Peaks_position = Peaks_data(ROI = ROI)
#%% Function to plot baseline
def Plot_baseline_fitting():
    fig, ax = plt.subplots(figsize = (10, 6))
    sns.set(style = 'ticks')
    ax.plot(Time_axis,
            Current_axis,
            label = 'Raw data',
            lw = 1.0,
            color = GrandBudapest5_5.hex_colors[0])
    ax.plot(Time_axis,
            baseline,
            label = 'SNIP fitting',
            lw = 2.0,
            alpha = 0.5,
            color = GrandBudapest5_5.hex_colors[1],
            linestyle = '--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim([0, int(math.ceil(np.max(Time_axis)))])
    ax.set_xlabel('Time (s)', fontsize = 18)
    ax.set_ylabel('Current (pA)', fontsize = 18)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.legend(fontsize = 16)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plt.tight_layout()
    plt.show()
    return
#%% Function to plot corrected baseline (centre at 0)
def Plot_baseline_corrected(Std_range_low = 1,
                            Std_range_high = 7):

    fig, ax = plt.subplots(figsize = (10, 6))
    sns.set(style = 'ticks')
    ax.plot(Time_axis,
            baseline_corrected,
            lw = 1.0,
            color = Mendl_4.hex_colors[3],
            )
    # Standard deviation values
    std_values = np.arange(Std_range_low, Std_range_high+1)
    # Loop over the standard deviation values and create axhline plots
    for std in std_values:
        ax.axhline(
            y = np.std(baseline_corrected) * -std, 
            linestyle = '--', 
            label = f'{std}σ', 
            lw = 2,
            color = Vivid_10.hex_colors[std-1], alpha=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim([0, int(math.ceil(np.max(Time_axis)))])
    ax.set_xlabel('Time (s)', fontsize = 18)
    ax.set_ylabel('Current (pA)', fontsize = 18)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.legend(fontsize = 16, bbox_to_anchor=(1.0, 0.5), loc = 'center left')
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plt.tight_layout()
    plt.show()
    return
#%% Function to plot the adjusted baselin (at 0) for the ROI version
def Plot_baseline_corrected_ROI(Std_range_low = 1,
                                Std_range_high = 7):
    df = pd.DataFrame({'Time_axis': Time_axis, 
                       'Current_baseline_adjusted': baseline_corrected})
    df_rows_selected = df[(df['Time_axis'] >= 4) & (df['Time_axis'] <= 9)]
    fig, ax = plt.subplots(figsize = (10, 6))
    sns.set(style = 'ticks')
    ax.plot(df_rows_selected['Time_axis'],
            df_rows_selected['Current_baseline_adjusted'],
            lw = 1.0,
            color = Mendl_4.hex_colors[3],
            )
    # Standard deviation values
    std_values = np.arange(Std_range_low, Std_range_high+1)
    # Loop over the standard deviation values and create axhline plots
    for std in std_values:
        ax.axhline(
            y = np.std(baseline_corrected) * -std, 
            linestyle = '--', 
            label = f'{std}σ',
            lw = 2,
            color = Vivid_10.hex_colors[std-1], alpha=0.5)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim([int(np.min(df_rows_selected['Time_axis'])), 
                 int(np.max(df_rows_selected['Time_axis']))
                 ])
    ax.set_xlabel('Time (s)', fontsize = 18)
    ax.set_ylabel('Current (pA)', fontsize = 18)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.legend(fontsize = 16, bbox_to_anchor=(1.0, 0.5), loc = 'center left')
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plt.tight_layout()
    plt.show()
    return
#%% Plot baseline
def Plot_baseline(ROI, Std_range_low=1, Std_range_high=7):
    if ROI == 'Yes':
        df = pd.DataFrame({'Time_axis': Time_axis, 
                           'Current_baseline_adjusted': baseline_corrected})
        df_rows_selected = df[(df['Time_axis'] >= 4) & (df['Time_axis'] <= 9)]
        title = 'Baseline Corrected (ROI)'
    elif ROI == 'No':
        df = pd.DataFrame({'Time_axis': df_baseline_corrected.iloc[:, 0], 
                           'Current_baseline_adjusted': df_baseline_corrected.iloc[:, 1]})
        df_rows_selected = df.copy()
        title = 'Baseline Corrected'
    else:
        print("Invalid value for ROI. Please choose 'Yes' or 'No'.")
        return
    
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.set(style='ticks')
    
    ax.plot(df_rows_selected['Time_axis'], df_rows_selected['Current_baseline_adjusted'], lw=1.0, color=Mendl_4.hex_colors[3])
    
    std_values = np.arange(Std_range_low, Std_range_high+1)
    for std in std_values:
        ax.axhline(y=np.std(df['Current_baseline_adjusted']) * -std, linestyle='--', label=f'{std}σ', lw=2,
                   color=Vivid_10.hex_colors[std-1], alpha=0.5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim([int(np.min(df_rows_selected['Time_axis'])), int(np.max(df_rows_selected['Time_axis']))])
    ax.set_xlabel('Time (s)', fontsize=18)
    ax.set_ylabel('Current (pA)', fontsize=18)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.legend(fontsize=16, bbox_to_anchor=(1.0, 0.5), loc='center left')
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plt.subplots_adjust(left=0.1, right=0.88, top=0.9, bottom=0.12)
    ax.set_title(title, fontsize=20, loc='left')
    plt.tight_layout()
    plt.show()
    return


#%% Plot onlye the traces preview to see the differnces between raw and adjusted
def Plot_Trace_preview():
    sns.set_theme(style = 'ticks')
    mpl.rcParams['axes.spines.left'] = True
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['axes.spines.bottom'] = True
    
    fig, axes = plt.subplots(2, 1, figsize = (12, 6), sharex = True)
    
    ax1 = axes[0]
    ax1.set_title('Raw trace', fontsize = 16)
    ax1.plot(Time_axis, 
             Current_axis, 
             color = GrandBudapest5_5.hex_colors[0], 
             linewidth=0.5)
    
    ax1.set_xlim([0, int(math.ceil(np.max(Time_axis)))])
    ax1.set_ylabel('Current (pA)', fontsize = 18)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax1.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    
    ax2 = axes[1]
    ax2.set_title('Baseline adjusted', fontsize = 16)
    ax2.plot(Time_axis,
             baseline_corrected,
             color = Mendl_4.hex_colors[3],
             linewidth=0.5)
    ax2.set_ylabel('Current (pA)', fontsize = 18)
    ax2.set_xlabel('Time (s)', fontsize = 18)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plt.tight_layout()
    plt.show()
    return

#%% Function to plot a trace with peaks identified
def Plot_peaks_pick(ROI):
    if ROI in ['Yes', 'No']:
        if ROI == 'No':
            df = df_baseline_corrected.copy()
        else:
            df = df_baseline_corrected_ROI.copy()
        Time = df.iloc[:, 0]
        peaks_signal = df.iloc[:, 1]
        peaks = df_Peak_search_result.iloc[:, 1]
        event_number = df_Peak_search_result.iloc[:, 0]
    else:
        print('Either Yes or No for the trace ROI')

    fig, ax = plt.subplots(figsize = (10, 6))
    sns.set(style = 'ticks')
    ax.plot(Time,
            peaks_signal,
            lw = 1.0,
            color = Mendl_4.hex_colors[3],
            )
    ax.plot(peaks/100000, peaks_signal[peaks], marker = 'o', markersize = 6, color = Mendl_4.hex_colors[0], linestyle = '')
    plt.show()
    for idx, (peak, event) in enumerate(zip(peaks, event_number)):
        ax.annotate(event, (peak/100000, peaks_signal[peak]))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim([int(np.min(Time)), 
                 int(math.ceil(np.max(Time)))
                 ])
    ax.set_xlabel('Time (s)', fontsize = 18)
    ax.set_ylabel('Current (pA)', fontsize = 18)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plt.tight_layout()
    plt.show()
    return

#%% Function to plot 20 random peaks
def Plot_Random_peaks_plot(Seed_range = 100,
                           Peak_numbers = 1,
                           Transparency = 1,
                           Color = Mendl_4.hex_colors[0]):
    # Generate a list of random numbers as seeds
    random_seeds = random.sample(range(Seed_range), 20)
    # Create a 5x5 grid of subplots
    fig, axes = plt.subplots(nrows=4, ncols=5, figsize=(16, 8))
    # Iterate over each subplot and plot the overlay
    for index, ax in enumerate(axes.flatten()):
        seed = random_seeds[index]
        random.seed(seed)
        selected_keys = random.sample(list(dfs_Peaks_position.keys()), Peak_numbers)
        # Iterate over the selected keys and plot the corresponding dataframes
        for key in selected_keys:
            df = dfs_Peaks_position[key]
            ax.plot(df['Time']*1000, df['Current'], color = Color, alpha=Transparency)
        #plot adjustment
        ax.set_xlim([0, 8])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        plt.tight_layout()
    
    plt.show()
    return

#%% Plot all translocation events as overlay
def Plot_overlay_all(Color = Mendl_4.hex_colors[0]):
    fig, ax = plt.subplots(figsize = (7,7))
    for i, (peak_number, peak_data) in enumerate(dfs_Peaks_position.items()):
        ax.plot(peak_data['Current'], 
                color = Color,
                alpha=0.1)
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        ax.set_xlabel('Time (µs)', fontsize = 18)
        ax.set_ylabel('Current (pA)', fontsize = 18)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        plt.tight_layout()
    plt.show()
    return

#%% This is to screen peaks overlay that best represent your dataset.
def Plot_Seed_screening(Seed_range = 100,
                        Peak_numbers = 20,
                        Transparency = 0.2,
                        Color = Mendl_4.hex_colors[0]):
    # Generate a list of random numbers as seeds
    random_seeds = random.sample(range(Seed_range), 50)
    # Create a 5x5 grid of subplots
    fig, axes = plt.subplots(nrows=5, ncols=10, figsize=(16, 8))
    # Iterate over each subplot and plot the overlay
    for index, ax in enumerate(axes.flatten()):
        seed = random_seeds[index]
        random.seed(seed)
        selected_keys = random.sample(list(dfs_Peaks_position.keys()), Peak_numbers)
        # Iterate over the selected keys and plot the corresponding dataframes
        for key in selected_keys:
            df = dfs_Peaks_position[key]
            ax.plot(df['Current'], color = Color, alpha=Transparency)
        # Remove labels and tick labels
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        # Show border
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        # Set plot title
        ax.set_title(f"Seed {seed}")
        plt.tight_layout()
    plt.show()
    return
#%% after identify a specific seed_number, the overlay can be plotted specifically
def Plot_Peaks_overlay(seed_number,
                       color = Mendl_4.hex_colors[0],
                       Transparency = 0.2):
    random.seed(seed_number)
    selected_keys = random.sample(list(dfs_Peaks_position.keys()), 20)
    # Create a new figure
    fig, ax = plt.subplots(figsize = (7,7))
    sns.set_theme(style = 'ticks')
    # Iterate over the selected keys and plot the corresponding dataframes
    for key in selected_keys:
        df = dfs_Peaks_position[key]
        ax.plot(df['Time']*1000, df['Current'], color = color, alpha = Transparency)

    # Show border
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.set_xlabel('Time (ms)', fontsize = 18)
    ax.set_ylabel('Current (pA)', fontsize = 18)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.set_xlim(0, 10)
    ax.set_ylim(-600, 100)
    plt.tight_layout()
    plt.show()
    return

#%% export the translocation data
excel_export(Path = folder_path, Name = File_name)

#%%
"""Below for lookin at how the fitting goes"""
Plot_baseline_fitting()
save_to_path_svg(Path = folder_path, Name = File_name + '_01_baseline_fitting')
save_to_path_png(Path = folder_path, Name = File_name + '_01_baseline_fitting')
plt.close('all')
#%%
"""Below for looking at how the baseline adjusted current looks"""
Plot_baseline(ROI = ROI)
save_to_path_svg(Path = folder_path, Name = File_name + '_02_adjusted_baseline')
save_to_path_png(Path = folder_path, Name = File_name + '_02_adjusted_baseline')
plt.close('all')
#%%
"""Below for looking at fitting and baseline adjusted current"""
Plot_Trace_preview()
save_to_path_svg(Path = folder_path, Name = File_name + '_03_baseline_and_adjusted')
save_to_path_png(Path = folder_path, Name = File_name + '_03_baseline_and_adjusted')
plt.close('all')
#%%
"""Below for peaks calling, what is consider as a peak"""
Plot_peaks_pick(ROI = ROI)
save_to_path_svg(Path = folder_path, Name = File_name + '_04_events_calling')
save_to_path_png(Path = folder_path, Name = File_name + '_04_events_calling')
plt.close('all')
#%%
"""Below for plotting 20 single peaks"""
Plot_Random_peaks_plot(Color='#DE8DB9')
save_to_path_svg(Path = folder_path, Name = File_name + '_05_20_random_events')
save_to_path_png(Path = folder_path, Name = File_name + '_05_20_random_events')
plt.close('all')
#%%
"""Below to overlay all the peaks"""
Plot_overlay_all()
save_to_path_svg(Path = folder_path, Name = File_name + '_06_all_events_overlay')
save_to_path_png(Path = folder_path, Name = File_name + '_06_all_events_overlay')
plt.close('all')
#%%
"""Below for screening the random seeds"""
Plot_Seed_screening(Color=Vivid_10.hex_colors[3])
save_to_path_svg(Path = folder_path, Name = File_name + '_07_seed_screening_for_overlay_plot')
save_to_path_png(Path = folder_path, Name = File_name + '_07_seed_screening_for_overlay_plot')
plt.close('all')
#%%
"""Below for plotting overlay, when you know a good seed number use the line below"""
Plot_Peaks_overlay(color = Vivid_10.hex_colors[3], seed_number = 26)
save_to_path_svg(Path = folder_path, Name = File_name + '_08_20_events_overlay')
save_to_path_png(Path = folder_path, Name = File_name + '_08_20_events_overlay')
plt.close('all')