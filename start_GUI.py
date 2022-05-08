import pandas as pd
import numpy as np
import tkinter as tk
from matplotlib import pyplot as plt    

# GUI---------
root = tk.Tk()
root.title('Signal Analysis')
# root.resizable(False,False)
root.geometry('1300x700')

# data_frame
    
data_frame = tk.LabelFrame(root, text='Import Data', padx=10, pady=10)
data_frame.grid(row=0, column=0, columnspan=2, padx=10, pady=10, sticky='nw')


def meas_drop(value):
    
    meas_list = list(df[comp_var.get()])
    
    global meas_var
    meas_var = tk.StringVar()
    meas_var.set('measurement')
    meas_drop = tk.OptionMenu(data_frame, meas_var, *meas_list)
    meas_drop.grid(row=2, column=0, padx=5, pady=5, sticky='w')

def readcsv():
    global df, t
    root.csv_name = tk.filedialog.askopenfilename(
        initialdir='/', title='Select A File',
        filetypes=(('csv files','*.csv'),('All files', '*.*')))
    df = pd.read_csv(root.csv_name, header=None)
    df.columns = [df.iloc[0], df.iloc[1]]  # set first and second rows as headers
    df = df[2:].reset_index(drop=True)
    # df.set_index('All calculations',inplace=True)  # set the time column as index
    t = df.iloc[:,-1]
    comp_list,w = zip(*list(df))
    comp_list = list(dict.fromkeys(comp_list))  # remove duplicates
    comp_list.pop()  # delete last column (time column)
    
    global comp_var    
    comp_var = tk.StringVar()
    comp_var.set('component')
    comp_drop = tk.OptionMenu(data_frame, comp_var, *comp_list,
                              command=meas_drop)
    comp_drop.grid(row=1, column=0, padx=5, pady=5, sticky='w')
    
    global info_label_plot
    info_label.destroy()
    info_label_plot = tk.Label(info_frame, 
                          text='After choosing component and measurement, click "plot"')
    info_label_plot.pack()
 
    
import_data = tk.Button(data_frame, text='Import data',command=readcsv)
import_data.grid(row=0, column=0, padx=5, pady=5, sticky='w')


import data_visualisation

plot_data = tk.Button(data_frame, text='Plot',
                      command=lambda: data_visualisation.plot_comp(
                          figure,plot_subframe,ax,df,t,comp_var,meas_var,
                          info_label_plot,info_frame))
plot_data.grid(row=3, column=0, padx=5, pady=5, sticky='w')


# plot_frame
plot_frame = tk.LabelFrame(root, text='Plot', padx=10, pady=10)
plot_frame.grid(row=0, column=1, rowspan=2, padx=10, pady=10, sticky='nw')

plot_subframe = tk.LabelFrame(plot_frame, text='', padx=1, pady=1)
plot_subframe.grid(row=1, column=0, padx=5, pady=5, sticky='nw')

figure = plt.Figure(figsize=(8,5))
ax = figure.add_subplot(111)


clear = tk.Button(plot_frame, text='Clear',
                  command=lambda: data_visualisation.clear_plot(plot_subframe,figure,ax))
clear.grid(row=0, column=0, padx=5, pady=5, sticky='nw')


# prony_analysis_frame
prony_analysis_frame = tk.LabelFrame(root, text='Prony Analysis', padx=10, pady=10)
prony_analysis_frame.grid(row=1, column=0, padx=10, pady=10, sticky='nw')

## time window sliding
tk.Label(prony_analysis_frame, text='Time windows sliding:').grid(
    row=0, column=0, sticky='w')

tk.Label(prony_analysis_frame, text='Start time [ms]:').grid(
    row=1, column=0, sticky='w')
start_time = tk.Entry(prony_analysis_frame)
start_time.grid(row=1, column=1, sticky='w')


tk.Label(prony_analysis_frame, text='Interval [ms]:').grid(
    row=2, column=0, sticky='w') 
interval = tk.Entry(prony_analysis_frame)
interval.grid(row=2, column=1, sticky='w')

tk.Label(prony_analysis_frame, text='Step [ms]:').grid(
    row=3, column=0, sticky='w') 
step = tk.Entry(prony_analysis_frame)
step.grid(row=3, column=1, sticky='w')

tk.Label(prony_analysis_frame, text='Stop time [ms]:').grid(
    row=4, column=0, sticky='w') 
stop_time = tk.Entry(prony_analysis_frame)
stop_time.grid(row=4, column=1, sticky='w')

## modal order
tk.Label(prony_analysis_frame, text='Modal order (n):').grid(
    row=5, column=0, sticky='w')

tk.Label(prony_analysis_frame, text='n_start:').grid(
    row=6, column=0, sticky='w') 
n_start = tk.Entry(prony_analysis_frame)
n_start.grid(row=6, column=1, sticky='w')

tk.Label(prony_analysis_frame, text='n_stop:').grid(
    row=7, column=0, sticky='w') 
n_stop = tk.Entry(prony_analysis_frame)
n_stop.grid(row=7, column=1, sticky='w')

tk.Label(prony_analysis_frame, text='n_step:').grid(
    row=8, column=0, sticky='w') 
n_step = tk.Entry(prony_analysis_frame)
n_step.grid(row=8, column=1, sticky='w')

tk.Label(prony_analysis_frame, text='n_break_criteria:').grid(
    row=9, column=0, sticky='w') 
n_break_criteria = tk.Entry(prony_analysis_frame)
n_break_criteria.grid(row=9, column=1, sticky='w')

## hard limit
tk.Label(prony_analysis_frame, text='Hard limits:').grid(
    row=10, column=0, sticky='w')

tk.Label(prony_analysis_frame, text='amp_limit_ratio:').grid(
    row=11, column=0, sticky='w') 
amp_limit_ratio = tk.Entry(prony_analysis_frame)
amp_limit_ratio.grid(row=11, column=1, sticky='w')

tk.Label(prony_analysis_frame, text='freq_limit:').grid(
    row=12, column=0, sticky='w') 
freq_limit = tk.Entry(prony_analysis_frame)
freq_limit.grid(row=12, column=1, sticky='w')

tk.Label(prony_analysis_frame, text='dr_limit:').grid(
    row=13, column=0, sticky='w') 
dr_limit = tk.Entry(prony_analysis_frame)
dr_limit.grid(row=13, column=1, sticky='w')

## run prony button
from prony_analysis import run_prony

prony_plot = tk.Button(prony_analysis_frame, text='Run Prony Analysis',
command= lambda: run_prony(root,
    int(start_time.get()),int(interval.get()),int(step.get()),int(stop_time.get()),
    int(n_start.get()),int(n_stop.get()),int(n_step.get()),float(n_break_criteria.get()),
    float(amp_limit_ratio.get()),float(freq_limit.get()),float(dr_limit.get()),
    t,df,comp_var,meas_var)
)

prony_plot.grid(row=14, column=0, sticky='w')

## prony default parameters
from prony_parameters import default_parameters
tk.Button(prony_analysis_frame, text='Default Parameters', command= lambda: default_parameters(
    start_time,interval,step,stop_time,n_start,n_stop,n_step,n_break_criteria,
    amp_limit_ratio,freq_limit,dr_limit)).grid(row=14, column=1, padx=5, pady=5, sticky='w')

## prony save parameters
from prony_parameters import save_parameters
tk.Button(prony_analysis_frame, text='Save Parameters', command=lambda: save_parameters(
    start_time,interval,step,stop_time,n_start,n_stop,n_step,n_break_criteria,
    amp_limit_ratio,freq_limit,dr_limit)).grid(row=15, column=1, padx=5, pady=5, sticky='w')

## prony import parameters
from prony_parameters import import_parameters
tk.Button(prony_analysis_frame, text='Import Parameters', command=lambda: import_parameters(
    start_time,interval,step,stop_time,n_start,n_stop,n_step,n_break_criteria,
    amp_limit_ratio,freq_limit,dr_limit)).grid(row=16, column=1, padx=5, pady=5, sticky='w')

## prony compare with other data
from compare_gui import compare
tk.Button(prony_analysis_frame, text='Compare Results', command=lambda: compare(
    root)).grid(row=15, column=0, sticky='w') 

# basic_analysis_frame
bsc_frame = tk.LabelFrame(root, text='Basic data analysis', padx=10, pady=10)
bsc_frame.grid(row=0, column=2, rowspan=3, padx=10, pady=10, sticky='nw')

window_b = tk.Label(bsc_frame, text='start analyse transience part \nfrom around time point [s]: ')#set start time for transient analysis
window_b.grid(row=2, column=0, sticky='w')
s = tk.Entry(bsc_frame)
s.grid(row=3, column=0, sticky='w')
window_e = tk.Label(bsc_frame, text='end at around time point [s]: ')#set start time for transient analysis
window_e.grid(row=4, column=0, sticky='w')
e = tk.Entry(bsc_frame)
e.grid(row=5, column=0, sticky='w')

from basic_analysis import dyn_analysis,trans,osci
trans_btn = tk.Button(bsc_frame, text='Transience part analysis', command= lambda: trans(df,s,e,comp_var,meas_var))
trans_btn.grid(row=6, column=0, sticky='nw')
osci_btn = tk.Button(bsc_frame, text='Oscillation part analysis', command= lambda: osci(df,comp_var,meas_var))
osci_btn.grid(row=7, column=0, sticky='nw')
dyn_btn = tk.Button(bsc_frame, text='Dynamic performance indices', command= lambda: dyn_analysis(df,comp_var,meas_var,bsc_frame,ax))
dyn_btn.grid(row=0, column=0, sticky='nw')

# multiple_analysis_frame
from multiple_analysis import readcsv_mlt, mlt_comp_dyn, mlt_comp_trans, mlt_comp_osci
mlt_frame = tk.LabelFrame(bsc_frame, text='Multiple Components analysis', padx=10, pady=10)
mlt_frame.grid(row=8, column=0, rowspan=3, padx=5, pady=5, sticky='nw')

msg1 = tk.Label(mlt_frame, text='1. Choose 1 file to analyse 1 \nmeasurement in all its components: ', anchor='w')
msg1.grid(row=0, column=0, sticky='w')

import_data_s = tk.Button(mlt_frame, text='Import data',command= lambda: readcsv_mlt(mlt_frame))
import_data_s.grid(row=0, column=1, padx=5, pady=5, sticky='w')

msg2 = tk.Label(mlt_frame, text='OR Choose multiple files to analyse \n1 measurement in all distances: ', anchor='w')
msg2.grid(row=1, column=0, sticky='w')

from multiple_analysis_2 import readcsv_mlt_dist
import_data_m = tk.Button(mlt_frame, text='Import multiple files',command= lambda: readcsv_mlt_dist(mlt_frame))
import_data_m.grid(row=1, column=1, padx=5, pady=5, sticky='w')

# info_frame
info_frame = tk.LabelFrame(root, text='Info/Hint', padx=10, pady=10)
info_frame.grid(row=3, column=0, columnspan=2, padx=5, pady=5, sticky='w')

info_label = tk.Label(info_frame, text='Import data')
info_label.pack()


# end of GUI--------
root.mainloop()
