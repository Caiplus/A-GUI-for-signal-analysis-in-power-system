import numpy as np
import pandas as pd
import tkinter as tk

def clear_all(start_time,interval,step,stop_time,n_start,n_stop,n_step,
              n_break_criteria,amp_limit_ratio,freq_limit,dr_limit):
    start_time.delete(0, 'end')
    interval.delete(0, 'end')
    step.delete(0, 'end')
    stop_time.delete(0, 'end')
    
    n_start.delete(0, 'end')
    n_stop.delete(0, 'end')
    n_step.delete(0, 'end')
    n_break_criteria.delete(0, 'end')
    
    amp_limit_ratio.delete(0, 'end')
    freq_limit.delete(0, 'end')
    dr_limit.delete(0, 'end')

def import_parameters(start_time,interval,step,stop_time,n_start,n_stop,n_step,
                       n_break_criteria,amp_limit_ratio,freq_limit,dr_limit):
    csv_name = tk.filedialog.askopenfilename(
        initialdir='/', title='Select A File',
        filetypes=(('csv files','*.csv'),('All files', '*.*')))
    df = pd.read_csv(csv_name, header=None)
    parameters = df[0].values.tolist()
    
    clear_all(start_time,interval,step,stop_time,n_start,n_stop,n_step,
                       n_break_criteria,amp_limit_ratio,freq_limit,dr_limit)
    
    start_time.insert(0, parameters[0])
    interval.insert(0, parameters[1])
    step.insert(0, parameters[2])
    stop_time.insert(0, parameters[3])
    
    n_start.insert(0, parameters[4])
    n_stop.insert(0, parameters[5])
    n_step.insert(0, parameters[6])
    n_break_criteria.insert(0, parameters[7])
    
    amp_limit_ratio.insert(0, parameters[8])
    freq_limit.insert(0, parameters[9])
    dr_limit.insert(0, parameters[10])

def save_parameters(start_time,interval,step,stop_time,n_start,n_stop,n_step,
                       n_break_criteria,amp_limit_ratio,freq_limit,dr_limit):
    prony_parameters = np.array([start_time.get(),interval.get(),step.get(),
                                 stop_time.get(),n_start.get(),n_stop.get(),
                                 n_step.get(),n_break_criteria.get(),amp_limit_ratio.get(),
                                 freq_limit.get(),dr_limit.get()]).astype(float)
    print(prony_parameters)
    np.savetxt('prony_parameters.csv', prony_parameters, delimiter=',')

def default_parameters(start_time,interval,step,stop_time,n_start,n_stop,n_step,
                       n_break_criteria,amp_limit_ratio,freq_limit,dr_limit):

    clear_all(start_time,interval,step,stop_time,n_start,n_stop,n_step,
                       n_break_criteria,amp_limit_ratio,freq_limit,dr_limit)
    
    start_time.insert(0,51000)
    interval.insert(0,10000)
    step.insert(0,1000)
    stop_time.insert(0,70000)
    
    n_start.insert(0,100)
    n_stop.insert(0,301)
    n_step.insert(0,40)
    n_break_criteria.insert(0,0.00001)
    
    amp_limit_ratio.insert(0,0.1)
    freq_limit.insert(0,10)
    dr_limit.insert(0,1)