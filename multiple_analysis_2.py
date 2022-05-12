import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import *
from tkinter import ttk
from tkinter.filedialog import asksaveasfile
from matplotlib import pyplot as plt
from scipy import signal
import seaborn as sns; sns.set_theme()
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.colors import LogNorm, Normalize
import csv

def options(mlt_frame):#the options of analysis
    msg2 = tk.Label(mlt_frame, text='2. Choose the value to be analysed')
    msg2.grid(row=2, column=0, sticky='w')
    
    msg3 = tk.Label(mlt_frame, text='3. Choose type of analysis')
    msg3.grid(row=3, column=0, sticky='w')
    
    mlt_dyn_btn = tk.Button(mlt_frame, text='Dynamic performance indices', 
                            command= lambda: mlt_comp_dyn(mlt_frame))
    mlt_dyn_btn.grid(row=4, column=0, sticky='nw')
    
    window_b_mlt = tk.Label(mlt_frame, 
                            text='start analyse transience part \nfrom around time point [s]: ')#set start time for transient analysis
    window_b_mlt.grid(row=5, column=0, sticky='w')
    s_mlt = tk.Entry(mlt_frame)
    s_mlt.grid(row=6, column=0, sticky='w')
    window_e_mlt = tk.Label(mlt_frame, 
                            text='end at around time point [s]: ')#set start time for transient analysis
    window_e_mlt.grid(row=7, column=0, sticky='w')
    e_mlt = tk.Entry(mlt_frame)
    e_mlt.grid(row=8, column=0, sticky='w')
    
    mlt_trans_btn = tk.Button(mlt_frame, text='Transience part analysis', 
                              command= lambda: mlt_comp_trans(mlt_frame, s_mlt, e_mlt))
    mlt_trans_btn.grid(row=9, column=0, sticky='nw')
    mlt_osci_btn = tk.Button(mlt_frame, text='Oscillation part analysis', 
                             command= lambda: mlt_comp_osci(mlt_frame, s_mlt, e_mlt))
    mlt_osci_btn.grid(row=10, column=0, sticky='nw')

def readcsv_mlt_dist(mlt_frame):
    global comp_list, meas_list, file_list
    #import files
    mlt_frame.csv_names = tk.filedialog.askopenfilenames(
        initialdir='/', title='Select Files',
        filetypes=(('csv files','*.csv'),('All files', '*.*')))
    file_list = list(mlt_frame.csv_names)
    comp_list = []
    meas_list = []
    for i in range(len(file_list)):
        df = pd.read_csv(file_list[i], header=None)
        df.columns = [df.iloc[0], df.iloc[1]]  # set first and second rows as headers
        df_data = df[2:].reset_index(drop=True)  #set from row 2 the main range of data
        comps,meass = zip(*list(df_data))  #from headers get the list of components and measurements
        comp_list.append(comps)
        meas_list.append(meass)
    comp_list = list(*dict.fromkeys(comp_list))  # remove duplicates
    comp_list.pop()  # delete last column (time column)
    meas_list = list(*dict.fromkeys(meas_list))
    meas_list = list(dict.fromkeys(meas_list)) # remove duplicates
    meas_list.pop()  # delete last column (time column)

    global meas_var
    meas_var = tk.StringVar()
    meas_var.set('measurement')
    meas_drop = tk.OptionMenu(mlt_frame, meas_var, *meas_list, command= options(mlt_frame))
    meas_drop.grid(row=2, column=1, padx=5, pady=5, sticky='w')

def mlt_dist_dyn(mlt_frame):
    mlt_dyn_top = tk.Toplevel()
    mlt_dyn_top.title('Multiple dynamic performance indices')
    for i in range(len(file_list)):  #loop for different files
        file = file_list[i]  #the file which is analysed
        df = pd.read_csv(file, header=None)
        df_data = df[2:].reset_index(drop=True)  #set from row 2 the main range of data
        #get the keywords of 'distance' from filename
        d = file.find('dist')  #locate the keyword 'dist'
        file2 = file[d:]
        l = file2.find('_')
        if l != -1:  #for there is '_' after distance number
            dist = file[d:d+l]
        else:  #fpr there is no '_' after distance number
            dist = file[d:-4]
        
        for column in range(len(df.columns)):
            if df.iloc[1,column] == meas_var.get():
                comp = df.iloc[0,column]
                t = df_data.iloc[:,-1].to_numpy(dtype=float)
                y = df_data[comp][meas_var.get()].to_numpy(dtype=float)
               
                steady_state = y[-1]#set the last value as steady state
                max_value = np.max(y)#find max value
                t_max = t[y.argmax()]#find time point of max value
                min_value = np.min(y)#find min value
                t_min = t[y.argmin()]#find time point of min value
                yr = y[y.argmin()+1:]#select the part after reaching the min value
                t_rise = t[np.where(yr>steady_state)[0][0]+y.argmin()+1]#find rise time on selected part
                osr = (max_value-steady_state)/steady_state#calculate overshoot ratio
    