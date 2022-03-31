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

def readcsv_mlt(mlt_frame):
    global df, df_data, t, comp_list, meas_list
    mlt_frame.csv_name = tk.filedialog.askopenfilename(
        initialdir='/', title='Select A File',
        filetypes=(('csv files','*.csv'),('All files', '*.*')))
    df = pd.read_csv(mlt_frame.csv_name, header=None)
    df.columns = [df.iloc[0], df.iloc[1]]  # set first and second rows as headers
    df_data = df[2:].reset_index(drop=True)
    t = df_data.iloc[:,-1]
    comp_list,meas_list = zip(*list(df_data))
    comp_list = list(dict.fromkeys(comp_list))  # remove duplicates
    comp_list.pop()  # delete last column (time column)
    meas_list = list(dict.fromkeys(meas_list))  # remove duplicates
    meas_list.pop()  # delete last column (time column)

    global meas_var
    meas_var = tk.StringVar()
    meas_var.set('measurement')
    meas_drop = tk.OptionMenu(mlt_frame, meas_var, *meas_list)
    meas_drop.grid(row=2, column=0, padx=5, pady=5, sticky='w')

def mlt_comp_dyn(mlt_frame):
    mlt_dyn_top = tk.Toplevel()
    mlt_dyn_top.title('Multiple dynamic performance indices')
    
    dyn_list = []
    hm_t_list = []
    hm_y_list = []  
    
    for column in range(len(df.columns)):
        if df.iloc[1,column] == meas_var.get():
            comp = df.iloc[0,column]
        #for comp in comp_list:
            t = df_data.iloc[:,-1].to_numpy(dtype=float)
            y = df_data[comp][meas_var.get()].to_numpy(dtype=float)
           
            steady_state = y[-1]
            max_value = np.max(y)
            t_max = t[y.argmax()]
            min_value = np.min(y)
            t_min = t[y.argmin()]
            yr = y[y.argmin()+1:]
            t_rise = t[np.where(yr>steady_state)[0][0]+y.argmin()+1]
            
            #for heat_map
            hm_t = [comp] + [t_max, t_min, t_rise]
            hm_t_list.append(hm_t)
            hm_y = [comp] + [steady_state, max_value, min_value]
            hm_y_list.append(hm_y)
            
            #write data to csv files
            dyn_data = [comp] + [steady_state, max_value, min_value, t_max, t_min, t_rise]
            dyn_list.append(dyn_data)
    
    header_dyn = [f'{meas_var.get()}', 'steady_state', 'max_value', 'min_value', 
                  't_max (s)', 't_min (s)', 't_rise (s)']
    with open('Dynamic performance indices (multiple).csv', 'w', newline='') as dyn:
        writer = csv.writer(dyn)
        writer.writerow(header_dyn)
        writer.writerows(dyn_list)

    header_dyn_y = [f'{meas_var.get()}', 'steady_state', 'max_value', 'min_value']
    with open('dyn_y.csv', 'w', newline='') as dyn_y:
        writer = csv.writer(dyn_y)
        writer.writerow(header_dyn_y)
        writer.writerows(hm_y_list) 
        
    header_dyn_t = [f'{meas_var.get()}', 't_max (s)', 't_min (s)', 't_rise (s)']
    with open('dyn_t.csv', 'w', newline='') as dyn_t:
        writer = csv.writer(dyn_t)
        writer.writerow(header_dyn_t)
        writer.writerows(hm_t_list)

    #build heat_map    
    hm_t_frame = tk.LabelFrame(mlt_dyn_top, text='', padx=10, pady=10)
    hm_t_frame.grid(row=0, column=1, padx=10, pady=10, sticky='nw')
    hm_y_frame = tk.LabelFrame(mlt_dyn_top, text='', padx=10, pady=10)
    hm_y_frame.grid(row=0, column=0, padx=10, pady=10, sticky='nw')
    
    data_y = pd.read_csv('dyn_y.csv',index_col=0)
    data_t = pd.read_csv('dyn_t.csv',index_col=0)
    
    figure1, ax = plt.subplots(figsize=(8,6))
    ax = sns.heatmap(data_y, square=True)
    ax.set_title('Heatmap of values', fontsize=16)
    
    figure2, ax = plt.subplots(figsize=(8,6))
    ax = sns.heatmap(data_t, square=True)
    ax.set_title('Heatmap of Timepoints', fontsize=16)
    
    canvas1 = FigureCanvasTkAgg(figure1, hm_y_frame)  # A tk.DrawingArea.
    canvas1.draw()
    canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    canvas2 = FigureCanvasTkAgg(figure2, hm_t_frame)  # A tk.DrawingArea.
    canvas2.draw()
    canvas2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas1, hm_y_frame)
    toolbar.update()
    canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    toolbar = NavigationToolbar2Tk(canvas2, hm_t_frame)
    toolbar.update()
    canvas2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
def mlt_comp_trans(mlt_frame, s_mlt, e_mlt):
    mlt_trans_top = tk.Toplevel()
    mlt_trans_top.title('Multiple Transience Analysis')
    
    trs_list = []
    hm_t_list = []
    hm_y_list = []
    hm_s_list = []
    hm_a_list = []
    
    for column in range(len(df.columns)):
        if df.iloc[1,column] == meas_var.get():
            comp = df.iloc[0,column]
            t = df_data.iloc[:,-1].to_numpy(dtype=float)
            y = df_data[comp][meas_var.get()].to_numpy(dtype=float)
            ##window from entry
            ts = float(s_mlt.get())
            ids = int(ts*1000)
            te = float(e_mlt.get())
            ide = int(te*1000)
            t = t[ids:ide]
            y = y[ids:ide]
            idt = np.argwhere(np.diff(y)).flatten()[0]#index of the last point before transient change
            ##first point P1(td,yd) after transient change
            td = t[idt+1]
            yd = y[idt+1]
            ##last point P2(tu,yu) before stepping up again
            tu = t[idt+150]
            yu = y[idt+150]
            ##calculate slope
            k = (yu-yd)/(tu-td)
            #calculate area
            Ar = np.sum(y[idt+1:idt+151])*(tu-td)/150
            
            trs_data = [comp] + [td, tu, yd, yu, k, Ar]
            trs_list.append(trs_data)
            hm_t_list.append([comp, td, tu])
            hm_y_list.append([comp, yd, yu])
            hm_s_list.append([comp, k])
            hm_a_list.append([comp, Ar])
            
    header_trs = [f'{meas_var.get()}', 't1 (s)', 't2 (s)', 'value1', 'value2', 'Slope', 'Area']    
    with open('Transient parameters (multiple).csv', 'w', newline='') as trs:
        writer = csv.writer(trs)
        writer.writerow(header_trs)
        writer.writerows(trs_list)
        
    header_trs_t = [f'{meas_var.get()}', 't1 (s)', 't2 (s)']    
    with open('trs_t.csv', 'w', newline='') as trs_t:
        writer = csv.writer(trs_t)
        writer.writerow(header_trs_t)
        writer.writerows(hm_t_list)
        
    header_trs_y = [f'{meas_var.get()}', 'value1', 'value2']    
    with open('trs_y.csv', 'w', newline='') as trs_y:
        writer = csv.writer(trs_y)
        writer.writerow(header_trs_y)
        writer.writerows(hm_y_list)
            
    header_trs_s = [f'{meas_var.get()}', 'Slope']    
    with open('trs_s.csv', 'w', newline='') as trs_s:
        writer = csv.writer(trs_s)
        writer.writerow(header_trs_s)
        writer.writerows(hm_s_list)

    header_trs_a = [f'{meas_var.get()}', 'Area']    
    with open('trs_a.csv', 'w', newline='') as trs_a:
        writer = csv.writer(trs_a)
        writer.writerow(header_trs_a)
        writer.writerows(hm_a_list)
        
    #build heat_map    
    hm_t_frame = tk.LabelFrame(mlt_trans_top, text='', padx=10, pady=10)
    hm_t_frame.grid(row=0, column=0, padx=10, pady=10, sticky='nw')
        
    hm_list = ['t', 'y', 's', 'a']
    for i in hm_list:
        n = hm_list.index(i)
        data = pd.read_csv(f'trs_{i}.csv',index_col=0)
        frame = tk.LabelFrame(hm_t_frame)
        frame.grid(row=0, column=n)
        
        figure, ax = plt.subplots(figsize=(4,6))
        ax = sns.heatmap(data, square=True, center=0)
        
        canvas = FigureCanvasTkAgg(figure, frame)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


        toolbar = NavigationToolbar2Tk(canvas, frame)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)       
        
def mlt_comp_osci(mlt_frame, s_mlt, e_mlt):
    mlt_osci_top = tk.Toplevel()
    mlt_osci_top.title('Multiple Oscillation Analysis')
    
    osci_list = []
    for column in range(len(df.columns)):
        if df.iloc[1,column] == meas_var.get():
            comp = df.iloc[0,column]
            t = df_data.iloc[:,-1].to_numpy(dtype=float)
            y = df_data[comp][meas_var.get()].to_numpy(dtype=float)        
            #define where Oscillation part starts (tb), ends (te)
            nadir_ind = signal.argrelextrema(y,np.less)[0]
            tb = np.max(nadir_ind)
            te = len(y)
            #Oscillation part
            x2 = t[tb:]
            y2 = y[tb:]
            #to draw line y=last_value_g
            g = y[-1]
            g2 = np.array([g] * (te-tb))
            #finding intersections of data plot and y=g
            idx = np.argwhere(np.diff(np.sign(y2 - g2))).flatten()#idx: relative indexes of intersections
            
            Area_list = []
            for i in range(len(idx)-1):
                n = idx[i]
                m = idx[i+1]
                ya = np.sum(y2[n:m+1])*(x2[m]-x2[n])/(m-n+1)
                ga = g*(x2[m]-x2[n])
                Area = abs(ya-ga)
                Area_list.append(Area)                
            osci_data = [comp] + [*Area_list]
            osci_list.append(osci_data)
    name_list = []
    for n in range(50):
        name_list.append('Area'+str(n))        
    header_osci = [f'{meas_var.get()}'] + name_list
    with open('Oscillation parameters (multiple).csv', 'w', newline='') as osci:
        writer = csv.writer(osci)
        writer.writerow(header_osci)
        writer.writerows(osci_list)
        
    #build heat_map    
    hm_a_frame = tk.LabelFrame(mlt_osci_top, text='', padx=10, pady=10)
    hm_a_frame.grid(row=0, column=0, padx=10, pady=10, sticky='nw')
        

    hm_a = pd.read_csv('Oscillation parameters (multiple).csv',index_col=0)
    
    figure, ax = plt.subplots(figsize=(12,6))
    ax = sns.heatmap(hm_a, norm=LogNorm())
    
    canvas = FigureCanvasTkAgg(figure, hm_a_frame)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, hm_a_frame)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)          
            
            
            
            
            
            