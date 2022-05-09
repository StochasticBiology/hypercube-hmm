import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl

mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18

sns.set_style("darkgrid", {"axes.facecolor": ".9"})




def bubble_plot(n_traits, p_matrix, names = [], x_label = 'Ordering', y_label = 'Trait', color = 'blue', header = '', fc = ''):
    
    if header[:2] == "Si":
        p_matrix = np.flipud(np.array(p_matrix))
        #print(p_matrix)
    
    x = []
    y = []
    s = []
    for i in range(n_traits):
        for j in range(n_traits):
            x.append(i)
            y.append(j)
            s.append(p_matrix[j][i]*2000)
    

    if fc == '':
        plt.scatter(x, y, s, c=color)
    else:
        plt.scatter(x, y, s, facecolors='none', edgecolors = color)
    
    plt.grid()
    plt.title(header, fontsize=20)
    plt.xlabel('Ordering', fontsize=18)
    plt.ylabel('Trait', fontsize=18)
    plt.xticks(list(range(n_traits)), [i+1 for i in list(range(n_traits))])
    if names == []:
        plt.yticks(list(range(n_traits)), list(range(n_traits)))
    else:
        plt.yticks(list(range(n_traits)), names)


def txt2matrix(txt_file):
    matrix = np.loadtxt(txt_file)
    return matrix


#Simple case plotting:

n_samples = [1]
for i in n_samples:
    bubble_plot(5, txt2matrix(f"mean_simple{i}.txt"), names = [1,2,3,4,5], color = 'dodgerblue', header = f'Simple case {i*4} samples (HBW)')
    bubble_plot(5, txt2matrix(f"sd_simple{i}.txt"), names = [1,2,3,4,5], color = 'black', fc = 'none', header = f'Simple case {i*4} samples (HBW)')
    plt.savefig(f"simple_{i}_samples_hyperHMM.png", dpi=900)
    plt.show()
    bubble_plot(5, txt2matrix(f"Plot files/hypertraps_single_L5_{i}.txt"), color = 'blue', header = f'HyperTraPS {i*4} samples')
    plt.grid()
    plt.savefig(f"simple_{i}_samples_hypertraps.png", dpi=900)
    plt.show()


