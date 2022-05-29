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
                            command= lambda: mlt_dist_dyn(mlt_frame))
    mlt_dyn_btn.grid(row=4, column=0, sticky='nw')
    
    window_b_mlt = tk.Label(mlt_frame, 
                            text='start analyse transient part \nfrom around time point [s]: ')#set start time for transient analysis
    window_b_mlt.grid(row=5, column=0, sticky='w')
    s_mlt = tk.Entry(mlt_frame)
    s_mlt.grid(row=6, column=0, sticky='w')
    window_e_mlt = tk.Label(mlt_frame, 
                            text='end at around time point [s]: ')#set start time for transient analysis
    window_e_mlt.grid(row=7, column=0, sticky='w')
    e_mlt = tk.Entry(mlt_frame)
    e_mlt.grid(row=8, column=0, sticky='w')
    
    mlt_trans_btn = tk.Button(mlt_frame, text='Transient part analysis', 
                              command= lambda: mlt_dist_trans(mlt_frame, s_mlt, e_mlt))
    mlt_trans_btn.grid(row=9, column=0, sticky='nw')
    mlt_osci_btn = tk.Button(mlt_frame, text='Oscillation part analysis', 
                             command= lambda: mlt_dist_osci(mlt_frame, s_mlt, e_mlt))
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
    for file in file_list:
        df = pd.read_csv(file, header=None)
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
    dyn_dist_top = tk.Toplevel()
    dyn_dist_top.title('Multiple dynamic performance indices')
    
    #loop for different files (distance)
    dyn_dist_list=[]
    ss_dist_list=[]
    tr_dist_list=[]
    osr_dist_list=[]
    for n in range(len(file_list)):  #the file which is analysed  
        df = pd.read_csv(file_list[n], header=None)
        df.columns = [df.iloc[0], df.iloc[1]]  # set first and second rows as headers
        df_data = df[2:].reset_index(drop=True)  #set from row 2 the main range of data
        
        #get the keywords of 'distance' from filename
        d = file_list[n].find('dist')  #locate the keyword 'dist'
        file2 = file_list[n][d:]
        if d != -1:  #for there is 'distx' in filename
            l = file2.find('_')
            if l != -1:  #for there is '_' after distance number
                dist = file_list[n][d:d+l]
            else:  #for there is no '_' after distance number
                dist = file_list[n][d:-4]
        else:  #for there is no 'distx' in filename
            dist = f'file{n}'
            
        #loop for different bus
        bus_list=[]
        steady_list=[]
        max_list=[]
        tmax_list=[]
        min_list=[]
        tmin_list=[]
        trise_list=[]
        osr_list=[]
        for column in range(len(df.columns)):
            if df.iloc[1,column] == meas_var.get():
                comp = df.iloc[0,column]
                t = df_data.iloc[:,-1].to_numpy(dtype=float)
                y = df_data[comp][meas_var.get()].to_numpy(dtype=float)
               
                steady_state = y[-1]#set the last value as steady state
                max_val = np.max(y)#find max value
                t_max = t[y.argmax()]#find time point of max value
                min_val = np.min(y)#find min value
                t_min = t[y.argmin()]#find time point of min value
                yr = y[y.argmin()+1:]#select the part after reaching the min value
                t_rise = t[np.where(yr>steady_state)[0][0]+y.argmin()+1]-50#find rise time on selected part
                osr = (max_val-steady_state)/steady_state#calculate overshoot ratio
                
                bus_list.append(comp)
                steady_list.append(steady_state)
                max_list.append(max_val)
                tmax_list.append(t_max)
                min_list.append(min_val)
                tmin_list.append(t_min)
                trise_list.append(t_rise)
                osr_list.append(osr)
        
        dyn_dist_val = [dist]+[*steady_list, *max_list, *min_list, *tmax_list, *tmin_list, 
                           *trise_list, *osr_list]
        dyn_dist_list.append(dyn_dist_val)
        ss_dist = [dist]+[*steady_list]
        ss_dist_list.append(ss_dist)
        tr_dist = [dist]+[*trise_list]
        tr_dist_list.append(tr_dist)
        osr_dist = [dist]+[*osr_list]
        osr_dist_list.append(osr_dist)
        
    #write datas to csv file
    dyn_btn = tk.Button(dyn_dist_top, text="save as csv", command = lambda: dyn_SaveFile())
    dyn_btn.grid(row=1, column=0, padx=10, pady=10, sticky='nw')
    def dyn_SaveFile():
        data = [("csv file(*.csv)","*.csv"),('All tyes(*.*)', '*.*')]
        dyn_file = tk.filedialog.asksaveasfilename(filetypes = data, defaultextension = data, 
                                                  initialfile = 'Dynamic performance indices for multiple distances.csv')
        # file will have file name provided by user.
        # Now we can use this file name to save file.
        n = len(bus_list)
        header_dyn1=['']+['steady_state']*n+['max_value']*n+['min_value']*n+['t_max(s)']*n+['t_min(s)']*n+['t_rise(s)']*n+['Overshoot Ratio']*n
        header_dyn2=['']+[*bus_list]*7
    
        with open(dyn_file, 'w', newline='') as dyn:
            writer = csv.writer(dyn)
            writer.writerow(header_dyn1)
            writer.writerow(header_dyn2)
            writer.writerows(dyn_dist_list)
            
    #write datas to the belonging csv file (for building heatmaps)
    header_ss = ['steady_state']+[*bus_list]
    with open('dyn_ss_dist.csv', 'w', newline='') as hm_ss:
            writer = csv.writer(hm_ss)
            writer.writerow(header_ss)
            writer.writerows(ss_dist_list)
    header_tr = ['rise_time(s)']+[*bus_list]
    with open('dyn_tr_dist.csv', 'w', newline='') as hm_tr:
            writer = csv.writer(hm_tr)
            writer.writerow(header_tr)
            writer.writerows(tr_dist_list)
    header_osr = ['overshoot_ratio']+[*bus_list]
    with open('dyn_osr_dist.csv', 'w', newline='') as hm_osr:
            writer = csv.writer(hm_osr)
            writer.writerow(header_osr)
            writer.writerows(osr_dist_list)
            
    #build heat_map    
    hm_dyn_frame2 = tk.LabelFrame(dyn_dist_top, text='', padx=10, pady=10)
    hm_dyn_frame2.grid(row=0, column=0, padx=10, pady=10, sticky='nw')
        
    hm_list = ['ss', 'tr', 'osr']
    for i in hm_list:
        n = hm_list.index(i)
        data = pd.read_csv(f'dyn_{i}_dist.csv',index_col=0)
        frame = tk.LabelFrame(hm_dyn_frame2)
        frame.grid(row=0, column=n)
        
        figure, ax = plt.subplots(figsize=(6,6))
        ax = sns.heatmap(data, annot=True)
        
        canvas = FigureCanvasTkAgg(figure, frame)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, frame)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1) 
    
def mlt_dist_trans(mlt_frame, s_mlt, e_mlt):
    trs_dist_top = tk.Toplevel()
    trs_dist_top.title('Multiple Transient Analysis')

    #loop for different files (distance)
    trs_dist_list=[]
    k_dist_list=[]
    ar_dist_list=[]
    for n in range(len(file_list)):  #the file which is analysed  
        df = pd.read_csv(file_list[n], header=None)
        df.columns = [df.iloc[0], df.iloc[1]]  # set first and second rows as headers
        df_data = df[2:].reset_index(drop=True)  #set from row 2 the main range of data
        
        #get the keywords of 'distance' from filename
        d = file_list[n].find('dist')  #locate the keyword 'dist'
        file2 = file_list[n][d:]
        if d != -1:  #for there is 'distx' in filename
            l = file2.find('_')
            if l != -1:  #for there is '_' after distance number
                dist = file_list[n][d:d+l]
            else:  #for there is no '_' after distance number
                dist = file_list[n][d:-4]
        else:  #for there is no 'distx' in filename
            dist = f'file{n}'
            
        #loop for different bus
        bus_list=[]
        td_list=[]
        tu_list=[]
        yd_list=[]
        yu_list=[]
        k_list=[]
        ar_list=[]
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
                
                bus_list.append(comp)
                td_list.append(td)
                tu_list.append(tu)
                yd_list.append(yd)
                yu_list.append(yu)
                k_list.append(k)
                ar_list.append(Ar)
                
        trs_dist_val = [dist]+[*td_list, *tu_list, *yd_list, *yu_list, *k_list, 
                       *ar_list]
        trs_dist_list.append(trs_dist_val)
        k_dist = [dist]+[*k_list]
        k_dist_list.append(k_dist)
        ar_dist = [dist]+[*ar_list]
        ar_dist_list.append(ar_dist)
        
    #write datas to csv file
    trs_btn = tk.Button(trs_dist_top, text="save as csv", command = lambda: trs_SaveFile())
    trs_btn.grid(row=1, column=0, padx=10, pady=10, sticky='nw')
    def trs_SaveFile():
        data = [("csv file(*.csv)","*.csv"),('All tyes(*.*)', '*.*')]
        trs_file = tk.filedialog.asksaveasfilename(filetypes = data, defaultextension = data, 
                                                  initialfile = 'Transient parameters for multiple distances.csv')
        # file will have file name provided by user.
        # Now we can use this file name to save file.
        n = len(bus_list)
        header_trs1=['']+['t_d(s)']*n+['t_u(s)']*n+['y_d']*n+['y_u']*n+['Slope']*n+['Area']*n
        header_trs2=['']+[*bus_list]*6
    
        with open(trs_file, 'w', newline='') as trs:
            writer = csv.writer(trs)
            writer.writerow(header_trs1)
            writer.writerow(header_trs2)
            writer.writerows(trs_dist_list)
                
    #write datas to the belonging csv file (for building heatmaps)
    header_k = ['Slope']+[*bus_list]
    with open('trs_k_dist.csv', 'w', newline='') as trs_k:
            writer = csv.writer(trs_k)
            writer.writerow(header_k)
            writer.writerows(k_dist_list)
    header_ar = ['Area']+[*bus_list]
    with open('trs_ar_dist.csv', 'w', newline='') as trs_ar:
            writer = csv.writer(trs_ar)
            writer.writerow(header_ar)
            writer.writerows(ar_dist_list)
            
    #build heat_map    
    hm_trs_frame2 = tk.LabelFrame(trs_dist_top, text='', padx=10, pady=10)
    hm_trs_frame2.grid(row=0, column=0, padx=10, pady=10, sticky='nw')
        
    hm_list = ['k', 'ar']
    for i in hm_list:
        n = hm_list.index(i)
        data = pd.read_csv(f'trs_{i}_dist.csv',index_col=0)
        frame = tk.LabelFrame(hm_trs_frame2)
        frame.grid(row=0, column=n)
        
        figure, ax = plt.subplots(figsize=(6,6))
        ax = sns.heatmap(data, center=0,  annot=True)
        
        canvas = FigureCanvasTkAgg(figure, frame)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, frame)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
def mlt_dist_osci(mlt_frame, s_mlt, e_mlt):
    osci_dist_top = tk.Toplevel()
    osci_dist_top.title('Multiple Oscillation Analysis')

    #loop for different files (distance)
    osci_list=[]
    Area_len=[]
    for n in range(len(file_list)):  #the file which is analysed  
        df = pd.read_csv(file_list[n], header=None)
        df.columns = [df.iloc[0], df.iloc[1]]  # set first and second rows as headers
        df_data = df[2:].reset_index(drop=True)  #set from row 2 the main range of data
        
        #get the keywords of 'distance' from filename
        d = file_list[n].find('dist')  #locate the keyword 'dist'
        file2 = file_list[n][d:]
        if d != -1:  #for there is 'distx' in filename
            l = file2.find('_')
            if l != -1:  #for there is '_' after distance number
                dist = file_list[n][d:d+l]
            else:  #for there is no '_' after distance number
                dist = file_list[n][d:-4]
        else:  #for there is no 'distx' in filename
            dist = f'file{n}'
            
        #loop for different bus
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
                    Area_len.append(len(Area_list))
                ar_sum = np.sum(Area_list)
                osci_data = [dist+' '+comp] + [ar_sum] + [*Area_list]
                osci_list.append(osci_data)
    #write the header according to number of areas
    n_ar = np.max(Area_len)
    name_list = []
    for n in range(n_ar):
        name_list.append('Area'+str(n))        
    header_osci = [f'{meas_var.get()}'] + ['Sum'] + name_list
    #write data to csv file
    osci_btn = tk.Button(osci_dist_top, text="save as csv", 
                         command = lambda: osci_SaveFile())
    osci_btn.grid(row=1, column=0, padx=10, pady=10, sticky='nw')
    def osci_SaveFile():
        data = [("csv file(*.csv)","*.csv"),('All tyes(*.*)', '*.*')]
        osci_mlt_file = tk.filedialog.asksaveasfilename(filetypes = data, defaultextension = data, 
                                                  initialfile = 'Oscillation parameters for multiple distances.csv')
        # file will have file name provided by user.
        # Now we can use this file name to save file.    
        with open(osci_mlt_file, 'w', newline='') as osci:
            writer = csv.writer(osci)
            writer.writerow(header_osci)
            writer.writerows(osci_list)
    
    #write datas to the belonging csv file (for building heatmaps)
    with open('osci_mlt.csv', 'w', newline='') as osci:
        writer = csv.writer(osci)
        writer.writerow(header_osci)
        writer.writerows(osci_list)
        
    #build heat_map    
    hm_osci_frame = tk.LabelFrame(osci_dist_top, text='', padx=10, pady=10)
    hm_osci_frame.grid(row=0, column=0, padx=10, pady=10, sticky='nw')        

    hm_osci = pd.read_csv('osci_mlt.csv',index_col=0)
    
    figure, ax = plt.subplots(figsize=(12,6))
    ax = sns.heatmap(hm_osci, vmin=0, norm=LogNorm())
    
    canvas = FigureCanvasTkAgg(figure, hm_osci_frame)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, hm_osci_frame)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
                    
                
                
                
    