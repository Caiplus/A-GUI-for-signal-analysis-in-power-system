import tkinter as tk
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def save_fig(title, figure):
    figure.savefig(title + '.png')

def plot_prony_compare(ax1,canvas):
    filename = tk.filedialog.askopenfilename(
    initialdir='/', title='Select A File',
    filetypes=(('csv files','*.csv'),('All files', '*.*')))
    modes = pd.read_csv(filename, header=None)
    modes = modes[1:].to_numpy(dtype=float)
    
    label = []
    for i in filename[::-1]:
        if i == '/':
            break
        label.append(i)
    label = label[::-1][:-4]
    print(''.join(label))
    ax1.scatter(modes[:,1], modes[:,2], label=''.join(label))
    ax1.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
    ax1.set_xlabel('Damping', fontsize=12)
    ax1.set_ylabel('Frequency [Hz]', fontsize=12)
    canvas.draw()


def compare(root):
    comp_gui = tk.Toplevel(root)
    comp_gui.title('Comparison between data')
    
    import_data = tk.Button(comp_gui, text='Import data',
                            command=lambda:plot_prony_compare(ax1,canvas))
    import_data.grid(row=0, column=0, padx=5, pady=5, sticky='w')
    
    figure = plt.Figure(figsize=(10,8))
    figure.set_tight_layout(True)
    ax1 = figure.add_subplot(111)
    ax1.set_title('Comparison between data', fontsize=14)
    
      
    canvas = FigureCanvasTkAgg(figure, comp_gui)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().grid(row=1, column=0, columnspan=2, sticky='w')
    
    title = 'Comparison between data'
    tk.Button(comp_gui, text='save as figure',
                            command=lambda:save_fig(title, figure)).grid(
                                row=0, column=1, padx=5, pady=5, sticky='e')