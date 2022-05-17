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
    
    mlt_dyn_btn = tk.Button(mlt_frame, text='Dynamic performance indices', command= lambda: mlt_comp_dyn(mlt_frame))
    mlt_dyn_btn.grid(row=4, column=0, sticky='nw')
    
    window_b_mlt = tk.Label(mlt_frame, text='start analyse transience part \nfrom around time point [s]: ')#set start time for transient analysis
    window_b_mlt.grid(row=5, column=0, sticky='w')
    s_mlt = tk.Entry(mlt_frame)
    s_mlt.grid(row=6, column=0, sticky='w')
    window_e_mlt = tk.Label(mlt_frame, text='end at around time point [s]: ')#set start time for transient analysis
    window_e_mlt.grid(row=7, column=0, sticky='w')
    e_mlt = tk.Entry(mlt_frame)
    e_mlt.grid(row=8, column=0, sticky='w')
    
    mlt_trans_btn = tk.Button(mlt_frame, text='Transience part analysis', command= lambda: mlt_comp_trans(mlt_frame, s_mlt, e_mlt))
    mlt_trans_btn.grid(row=9, column=0, sticky='nw')
    mlt_osci_btn = tk.Button(mlt_frame, text='Oscillation part analysis', command= lambda: mlt_comp_osci(mlt_frame, s_mlt, e_mlt))
    mlt_osci_btn.grid(row=10, column=0, sticky='nw')

def readcsv_mlt(mlt_frame):
    global dist, df, df_data, comp_list, meas_list
    mlt_frame.csv_name = tk.filedialog.askopenfilename(
        initialdir='/', title='Select A File',
        filetypes=(('csv files','*.csv'),('All files', '*.*')))
    
    #get the keywords of 'distance' from filename
    d = mlt_frame.csv_name.find('dist')  #locate the keyword 'dist'
    file2 = mlt_frame.csv_name[d:]
    if d != -1:  #for there is 'distx' in filename
        l = file2.find('_')
        if l != -1:  #for there is '_' after distance number
            dist = mlt_frame.csv_name[d:d+l]
        else:  #for there is no '_' after distance number
            dist = mlt_frame.csv_name[d:-4]
    else:  #for there is no 'distx' in filename
        dist = '_'
        
    #dataframe    
    df = pd.read_csv(mlt_frame.csv_name, header=None)
    df.columns = [df.iloc[0], df.iloc[1]]  # set first and second rows as headers
    df_data = df[2:].reset_index(drop=True)
    comp_list,meas_list = zip(*list(df_data))
    comp_list = list(dict.fromkeys(comp_list))  # remove duplicates
    comp_list.pop()  # delete last column (time column)
    meas_list = list(dict.fromkeys(meas_list))  # remove duplicates
    meas_list.pop()  # delete last column (time column)

    global meas_var
    meas_var = tk.StringVar()
    meas_var.set('measurement')
    meas_drop = tk.OptionMenu(mlt_frame, meas_var, *meas_list, command= options(mlt_frame))
    meas_drop.grid(row=2, column=1, padx=5, pady=5, sticky='w')

def mlt_comp_dyn(mlt_frame):
    mlt_dyn_top = tk.Toplevel()
    mlt_dyn_top.title('Multiple dynamic performance indices')
    
    dyn_list = []
    hm_t_list = []
    hm_y_list = []
    hm_osr_list = []
    
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
            
            #for heat_map
            hm_t = [comp] + [t_max, t_min, t_rise]#pick out all time values to make comparision
            hm_t_list.append(hm_t)#put in list
            hm_y = [comp] + [steady_state, max_value, min_value]#pick out all goal values for comparision
            hm_y_list.append(hm_y)
            hm_osr = [comp] + [osr]#pick out overshoot ratio for comparision
            hm_osr_list.append(hm_osr)
            
            #write all data to csv file
            dyn_data = [comp] + [steady_state, max_value, min_value, t_max, t_min, t_rise, osr]
            dyn_list.append(dyn_data)
    
    #write datas to csv file
    dyn_btn = tk.Button(mlt_dyn_top, text="save as csv", command = lambda: dyn_SaveFile())
    dyn_btn.grid(row=1, column=0, padx=10, pady=10, sticky='nw')
    def dyn_SaveFile():
        data = [("csv file(*.csv)","*.csv"),('All tyes(*.*)', '*.*')]
        dyn_file = tk.filedialog.asksaveasfilename(filetypes = data, defaultextension = data, 
                                                  initialfile = f'Multiple Dynamic performance indices_{dist}.csv')
        # file will have file name provided by user.
        # Now we can use this file name to save file.
        header_dyn = [f'{meas_var.get()}', 'steady_state', 'max_value', 'min_value', 
                      't_max (s)', 't_min (s)', 't_rise (s)', 'Overshoot Ratio']
        with open(dyn_file, 'w', newline='') as dyn:
            writer = csv.writer(dyn)
            writer.writerow(header_dyn)
            writer.writerows(dyn_list)
            
   #write datas to the belonging csv file (for building heatmaps)
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

    header_dyn_osr = [f'{meas_var.get()}', 'Overshoot Ratio']
    with open('dyn_osr.csv', 'w', newline='') as dyn_osr:
        writer = csv.writer(dyn_osr)
        writer.writerow(header_dyn_osr)
        writer.writerows(hm_osr_list)
    
    #build heat_map    
    hm_dyn_frame = tk.LabelFrame(mlt_dyn_top, text='', padx=10, pady=10)
    hm_dyn_frame.grid(row=0, column=0, padx=10, pady=10, sticky='nw')
        
    hm_list = ['y', 'osr', 't']
    for i in hm_list:
        n = hm_list.index(i)
        data = pd.read_csv(f'dyn_{i}.csv',index_col=0)
        frame = tk.LabelFrame(hm_dyn_frame)
        frame.grid(row=0, column=n)
        
        figure, ax = plt.subplots(figsize=(6,6))
        ax = sns.heatmap(data, square=True,  annot=True)
        
        canvas = FigureCanvasTkAgg(figure, frame)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, frame)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1) 
    
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

    #write data to csv file
    trs_btn = tk.Button(mlt_trans_top, text="save as csv", command = lambda: trans_SaveFile())
    trs_btn.grid(row=1, column=0, padx=10, pady=10, sticky='nw')
    def trans_SaveFile():
        data = [("csv file(*.csv)","*.csv"),('All tyes(*.*)', '*.*')]
        trans_file = tk.filedialog.asksaveasfilename(filetypes = data, defaultextension = data, 
                                                  initialfile = f'Multiple Transient parameters_{dist}.csv')
        # file will have file name provided by user.
        # Now we can use this file name to save file.            
        header_trs = [f'{meas_var.get()}', 't1 (s)', 't2 (s)', 'value1', 'value2', 'Slope', 'Area']    
        with open(trans_file, 'w', newline='') as trs:
            writer = csv.writer(trs)
            writer.writerow(header_trs)
            writer.writerows(trs_list)
            
   #write datas to the belonging csv file (for building heatmaps)        
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
        ax = sns.heatmap(data, square=True, center=0,  annot=True)
        
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
            Area_len = []
            for i in range(len(idx)-1):
                n = idx[i]
                m = idx[i+1]
                ya = np.sum(y2[n:m+1])*(x2[m]-x2[n])/(m-n+1)
                ga = g*(x2[m]-x2[n])
                Area = abs(ya-ga)
                Area_list.append(Area)
                Area_len.append(len(Area_list))
            osci_data = [comp] + [*Area_list]
            osci_list.append(osci_data)
    #write the header according to number of areas
    name_list = []
    n_ar = np.max(Area_len)+1
    for n in range(n_ar):
        name_list.append('Area'+str(n))        
    header_osci = [f'{meas_var.get()}'] + name_list

    #write data to csv file
    osci_btn = tk.Button(mlt_osci_top, text="save as csv", 
                         command = lambda: osci_SaveFile())
    osci_btn.grid(row=1, column=0, padx=10, pady=10, sticky='nw')
    def osci_SaveFile():
        data = [("csv file(*.csv)","*.csv"),('All tyes(*.*)', '*.*')]
        osci_file = tk.filedialog.asksaveasfilename(filetypes = data, defaultextension = data, 
                                                  initialfile = f'Multiple Oscillation parameters_{dist}.csv')
        # file will have file name provided by user.
        # Now we can use this file name to save file.    
        with open(osci_file, 'w', newline='') as osci:
            writer = csv.writer(osci)
            writer.writerow(header_osci)
            writer.writerows(osci_list)
    
    #write datas to the belonging csv file (for building heatmaps)
    with open('osci_ar.csv', 'w', newline='') as osci:
        writer = csv.writer(osci)
        writer.writerow(header_osci)
        writer.writerows(osci_list)
        
    #build heat_map    
    hm_a_frame = tk.LabelFrame(mlt_osci_top, text='', padx=10, pady=10)
    hm_a_frame.grid(row=0, column=0, padx=10, pady=10, sticky='nw')        

    hm_a = pd.read_csv('osci_ar.csv',index_col=0)
    
    figure, ax = plt.subplots(figsize=(12,6))
    ax = sns.heatmap(hm_a, norm=LogNorm())
    
    canvas = FigureCanvasTkAgg(figure, hm_a_frame)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, hm_a_frame)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)          
            
            
            
            
            
            