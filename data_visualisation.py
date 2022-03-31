from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
import tkinter as tk


def plot_comp(figure,plot_subframe,ax,df,t,comp_var,meas_var):
    
    ax.plot(t.to_numpy(dtype=float), df[comp_var.get()][meas_var.get()].to_numpy(dtype=float),label=meas_var.get())
    ax.set_ylabel(meas_var.get(), fontsize=12)

    ax.set_title(comp_var.get() + meas_var.get(), fontsize=16)   
    ax.set_xlabel('time [$s$]', fontsize=12)
    ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper right')


    canvas = FigureCanvasTkAgg(figure, plot_subframe)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, plot_subframe)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    def on_key_press(event):
        print("you pressed {}".format(event.key))
        key_press_handler(event, canvas, toolbar)
    
    canvas.mpl_connect("key_press_event", on_key_press)

        
def clear_plot(plot_subframe):
    for widget in plot_subframe.winfo_children():
        widget.destroy()