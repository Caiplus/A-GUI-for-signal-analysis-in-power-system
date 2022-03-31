from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import numpy.polynomial.polynomial as poly
import cmath
import tkinter as tk
from tkinter import ttk
from scipy.linalg import toeplitz
from sklearn.cluster import DBSCAN

def save_fig(title, figure_p):
    figure_p.savefig(title + '.png')

def save_csv(title, modes):
    np.savetxt(f'{title}.csv', modes, delimiter=',', header='mode,damping_ratio,frequency')

class T_win:
    def __init__(self, t_i, y, t, interval,
                 n_start,n_stop,n_step,n_break_criteria,amp_limit_ratio,freq_limit,dr_limit,
                 scrollable_frame, column_index):
        
        self.t_win = t.iloc[t_i:t_i+interval].to_numpy(dtype=float)
        self.y_k = y.iloc[t_i:t_i+interval].to_numpy(dtype=float)
        
    # def prony_n(self, n_start,n_stop,n_step,n_break_criteria,amp_limit_ratio,freq_limit):
        
        for n in range(n_start,n_stop,n_step):
        
            # Apply prony analysis
            self.B_i, self.lamda_i, self.y_est = prony(self.t_win, self.y_k, n)
            
            # Mean squared error
            self.MSE = np.square(np.subtract(self.y_k,self.y_est)).mean().real
            
            # damping ratio
            self.dr = -self.lamda_i.real/abs(self.lamda_i)
            self.freq = abs(self.lamda_i.imag)/(2*cmath.pi)
            
            # mode shape
            self.amplitude = 2*abs(self.B_i)
            self.phase = np.angle(self.B_i)
            
            self.prony_shape = np.vstack((self.dr,self.freq,self.amplitude,self.phase,
                                          self.lamda_i,self.B_i)).transpose()
            
            # hard limit
            # self.amp_limit = np.amax(self.prony_shape[:,2]) *amp_limit_ratio  # relative to highest amplitude
            self.amp_limit = np.partition(self.prony_shape[:,2], -2)[-2] *amp_limit_ratio  # relative to second highest amplitude
            
            self.prony_shape_hl = np.delete(self.prony_shape,np.where(self.prony_shape[:,2]<self.amp_limit),axis=0)
            self.prony_shape_hl = np.delete(self.prony_shape_hl,np.where(self.prony_shape_hl[:,1]>freq_limit),axis=0)
            self.prony_shape_hl = np.delete(self.prony_shape_hl,np.where(self.prony_shape_hl[:,0]>dr_limit),axis=0)
            # self.prony_shape_hl = np.delete(self.prony_shape_hl,np.where(self.prony_shape_hl[:,0]<0),axis=0)

            
            self.n_hl = len(self.prony_shape_hl[:,2])
     
            # reconstruct the signal after hard limit:
            self.y_hl = [sum(self.prony_shape_hl[n_i,5]*cmath.exp(self.prony_shape_hl[n_i,4]*(t_k-self.t_win[0])) 
                        for n_i in range(self.n_hl)).real
                      for t_k in self.t_win]
            
            # Mean squared error after hard limit
            self.MSE_hl = np.square(np.subtract(self.y_k,self.y_hl)).mean().real
                    
            # prony_shape = np.vstack([['DR','Freq','Amplitude','Phase'], prony_shape])
            
            if self.MSE <= n_break_criteria:
                break
    
        
        # # plot prony fitting for the time window
        # plt.figure(column_index)
        # plt.plot(self.t_win, self.y_k, label=f'original data, {t_i/1000}')
        # plt.plot(self.t_win, self.y_est, linestyle='dashed', label=f'n={n}\nMSE={self.MSE}')
        # plt.plot(self.t_win, self.y_hl, linestyle='dotted', 
        #           label=f'After hard limit\nn={self.n_hl}\nMSE={self.MSE_hl}')
        # plt.xlabel('t', fontsize=12)
        # plt.ylabel('y(t)', fontsize=12)    
        # plt.title(f'Time window: {t_i} - {t_i+interval}', fontsize=16) 
        # plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
    
    # def t_win_plot(self, scrollable_frame):
        self.figure_p = plt.Figure(figsize=(5,3))
        self.ax = self.figure_p.add_subplot(111) 
        self.ax.plot(self.t_win, self.y_k, label='original data')
        self.ax.plot(self.t_win, self.y_est, linestyle='dashed', label=f'n={n}\nMSE={self.MSE}')
        self.ax.plot(self.t_win, self.y_hl, linestyle='dotted', 
                  label=f'After hard limit\nn={self.n_hl}\nMSE={self.MSE_hl}')
        
        self.ax.set_title(f'Time window: {t_i} - {t_i+interval}', fontsize=14)
        self.ax.set_xlabel('t', fontsize=12)
        self.ax.set_ylabel('y(t)', fontsize=12)
        self.ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper right')
        
        self.canvas1 = FigureCanvasTkAgg(self.figure_p, scrollable_frame)  # A tk.DrawingArea.
        self.canvas1.draw()
        self.canvas1.get_tk_widget().grid(row=0, column=column_index, sticky='nw')
        
        self.var1 = tk.IntVar()
        self.check = tk.Checkbutton(scrollable_frame, text='Well fitting',
                                    variable=self.var1, onvalue=1)
        self.check.grid(row=1, column=column_index, sticky='nw')


def prony(t, y_k, n):
 	# Step 1: to assemble the data into a Toeplitz matrix, Tmat
    # Tmat is (k-n)*n and bmat is (k-n)*1
    
    k = len(t)
    bmat = y_k[n:k]
    Tmat = toeplitz(y_k[n-1:k-1] , y_k[0:n][::-1])
	
    # Step 2: to solve polynomial coefficients a using least squares solution
    a = np.linalg.lstsq(Tmat, bmat, rcond=None)[0]

	# Step 3: to solve zi (roots of the polynomial of(8))
	# first, form the polynomial coefficients c
    c = np.zeros((n+1), dtype = 'complex_')
    c[n] = 1.
    for i in range(1,n+1):
        c[n-i] = -a[i-1]

    zi = poly.polyroots(c)
    
    # Step 4: find lamda_i and B_i to reconstruct y_est
    lamda_i = np.log(zi)/(t[1]-t[0])  # eigenvalues

	# Vandermonde matrix
    Vmat = np.zeros((k, n), dtype = 'complex_')
    for irow in range(k):
        Vmat[irow, :] = zi**irow
        
    B_i = np.linalg.lstsq(Vmat, y_k, rcond=None)[0]  #residues
    
    # reconstruct the signal
    y_est = [sum(B_i[n_i]*cmath.exp(lamda_i[n_i]*(t_k-t[0])) 
                      for n_i in range(n)).real
              for t_k in t]
    
    return B_i, lamda_i, y_est


def run_cluster(eps, min_samples, ax2, canvas, prony_gui, checkbox_list, var_chooseall, figure_p, root):
    cluster = np.empty((0,2), float)
    
    for win_t_i in checkbox_list:
        if globals()[win_t_i].var1.get() == 1 and var_chooseall.get() == 0:
            cluster = np.append(cluster, globals()[win_t_i].prony_shape_hl[:,0:2].real,axis=0)
        
        if var_chooseall.get() == 1:
            cluster = np.append(cluster, globals()[win_t_i].prony_shape_hl[:,0:2].real,axis=0)

    # normalization of cluster data using z-score scaling
    cluster_offset = cluster - np.average(cluster,axis=0)
    cluster_nor = cluster_offset/np.sqrt(np.var(cluster_offset,axis=0))
    
    clustering = DBSCAN(eps=eps,min_samples=min_samples).fit(cluster_nor)
    labels = clustering.labels_
    analysis = np.hstack((cluster,np.asarray(labels).reshape((len(cluster), 1))))
    
    colors = list(map(lambda x: 'y' if x == -1 else 'r', labels))
    
    modes = np.empty((0,2), float)
    
    n_modes = lambda x: len(set(labels))-1 if -1 in labels else len(set(labels))
    len_modes = n_modes(labels)
    
    for labels_i in range(len_modes):
        mode_i = np.average(cluster[np.where(analysis[:,2]==labels_i)], axis=0)
        modes = np.append(modes,[mode_i], axis=0)
    
    # scatter plot for damping and freq distribution after clustering
    ax2.scatter(cluster[:,0], cluster[:,1], c=colors, marker="o", picker=True)
    ax2.scatter(modes[:,0], modes[:,1], marker='x', color='#000000')
    ax2.set_xlabel('Damping ratio', fontsize=12)
    ax2.set_ylabel('Frequency [Hz]', fontsize=12)
    for xy in modes:
        ax2.annotate(f'({np.round(xy[0],3)},{np.round(xy[1],3)})', xy=xy, 
                      textcoords='data')
    canvas.draw()
    
    
    mode_index = np.array([i for i in range(len_modes)]).reshape(len_modes,1)
    modes = np.hstack((mode_index, modes))
    
    for row in modes:
        print(' '.join(map(lambda x: "{:.3f}\t".format(x), row)))

    print(modes)

    # save results
    tk.Button(prony_gui, text='save as figure',command= lambda: save_fig(
        title, figure_p)).grid(row=2, column=3, sticky='w')
    
    tk.Button(prony_gui, text='save modes as csv', command=lambda: save_csv(
        title, modes)).grid(row=3, column=3, sticky='w')



def run_prony(root, start_time, interval, step, stop_time,
          n_start, n_stop, n_step, n_break_criteria,
          amp_limit_ratio, freq_limit, dr_limit,
          t, df, comp_var, meas_var):     
    
    y = df[comp_var.get()][meas_var.get()]
    global title
    title = comp_var.get() + meas_var.get()
    

    prony_gui = tk.Toplevel(root)
    prony_gui.title('Prony Analysis')
    
    # well-fitting time window frame
    time_win_frame = tk.LabelFrame(prony_gui, text='Choose the well-fitting time windows:',
                                    padx=10, pady=10, height=10, width=50)
    time_win_frame.grid(row=0, column=0, columnspan=5, padx=10, pady=10, sticky='nw')
    time_win_frame.grid_propagate(0)
    
    canvas_twin = tk.Canvas(time_win_frame)
    scrollbar = ttk.Scrollbar(time_win_frame, orient="horizontal", command=canvas_twin.xview)
    scrollable_frame = tk.LabelFrame(canvas_twin)
    canvas_twin.pack(side="top", fill="both", expand=True)
    scrollbar.pack(side="bottom", fill="x")
    
    scrollable_frame.bind(
        "<Configure>",
        lambda e: canvas_twin.configure(
            scrollregion=canvas_twin.bbox("all")
        )
    )
    
    canvas_twin.create_window((0, 0), window=scrollable_frame, anchor="nw")
    
    canvas_twin.configure(xscrollcommand=scrollbar.set)

    var_chooseall = tk.IntVar()
    tk.Checkbutton(time_win_frame, text='Choose all time windows',
                                variable=var_chooseall, onvalue=1).pack()


    
    figure_p = plt.Figure(figsize=(10,5))
    figure_p.set_tight_layout(True)
    ax1 = figure_p.add_subplot(121)       # freq-damping scatter plot
    ax1.set_title('Before Clustering', fontsize=14)
    ax2 = figure_p.add_subplot(122)       # after clustering scatter plot
    ax2.set_title('After Clustering', fontsize=14)
      
    canvas = FigureCanvasTkAgg(figure_p, prony_gui)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().grid(row=1, column=0, columnspan=4, sticky='w')


    column_index = 0
    checkbox_list = []
    for t_i in range(start_time,stop_time,step):
        
        checkbox_list.append('win_'+str(t_i))
        
        globals()['win_'+str(t_i)] = T_win(t_i, y, t, interval,
                 n_start,n_stop,n_step,n_break_criteria,amp_limit_ratio,freq_limit,dr_limit,
                 scrollable_frame, column_index)

        column_index += 1
        
        # scatter plot for damping and freq distribution
        ax1.scatter(globals()['win_'+str(t_i)].prony_shape_hl[:,0], globals()['win_'+str(t_i)].prony_shape_hl[:,1], label=t_i/1000)
        ax1.set_xlabel('Damping', fontsize=12)
        ax1.set_ylabel('Freq [Hz]', fontsize=12)
        ax1.legend(bbox_to_anchor=(1.0, 1.0), loc='upper right')
        
        canvas.draw()
    
    
    # parameters input for cluster method
    tk.Label(prony_gui, text='eps:').grid(row=2, column=0, sticky='w')
    eps = tk.Entry(prony_gui)
    eps.grid(row=2, column=1, sticky='w')
    
    tk.Label(prony_gui, text='min_samples:').grid(row=3, column=0, sticky='w')
    min_samples = tk.Entry(prony_gui)
    min_samples.grid(row=3, column=1, sticky='w')
    
    # run cluster method
    tk.Button(prony_gui, text='Run Cluster Method',command= lambda: run_cluster(
        float(eps.get()),int(min_samples.get()),ax2,canvas,prony_gui,checkbox_list,
        var_chooseall,figure_p,root)).grid(row=2, column=2, sticky='w')
    

