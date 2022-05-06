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

n_samples = [1, 2, 4, 10]
for i in n_samples:
    bubble_plot(5, txt2matrix(f"Plot files/mean_simple{i}.txt"), names = [1,2,3,4,5], color = 'dodgerblue', header = f'Simple case {i*4} samples (HBW)')
    bubble_plot(5, txt2matrix(f"Plot files/sd_simple{i}.txt"), names = [1,2,3,4,5], color = 'black', fc = 'none', header = f'Simple case {i*4} samples (HBW)')
    plt.savefig(f"simple_{i}_samples_hyperHMM.png", dpi=900)
    plt.show()
    bubble_plot(5, txt2matrix(f"Plot files/hypertraps_single_L5_{i}.txt"), color = 'blue', header = f'HyperTraPS {i*4} samples')
    plt.grid()
    plt.savefig(f"simple_{i}_samples_hypertraps.png", dpi=900)
    plt.show()





#Ovarian cancer plotting:

bubble_plot(7, txt2matrix("Plot files/mean_ovarian.txt"), names = ['Xp-', '1q+', '8p-', '4q-', '5q-', '3q+', '8q+'], color = 'dodgerblue', header = 'Ovarian cancer data (HBW)')
bubble_plot(7, txt2matrix("Plot files/sd_ovarian.txt"), names = ['Xp-', '1q+', '8p-', '4q-', '5q-', '3q+', '8q+'], color = 'black', fc = 'none', header = 'Ovarian cancer data (HBW)')
plt.savefig(f"ovarian_plot_hyperHMM.png", dpi=900)
plt.show()

bubble_plot(7, txt2matrix(f"Plot files/hypertraps_ovarian.txt"),  names = ['Xp-', '1q+', '8p-', '4q-', '5q-', '3q+', '8q+'], color = 'blue', header = f'Ovarian cancer data (HyperTraPS)')
plt.grid()
plt.savefig(f"ovarian_plot_hypertraps.png", dpi=900)
plt.show()




#TB-drug plotting:

bubble_plot(10, txt2matrix("Plot files/mean_tb_drug.txt"), names = ["AMI", "CAP", "MOX", "OFL", "PRO", "PZA", "EMB", "STR", "RIF", "INH"], color = 'dodgerblue', header = 'TB-drub data (HBW)')
bubble_plot(10, txt2matrix("Plot files/sd_tb_drug.txt"), names = ["AMI", "CAP", "MOX", "OFL", "PRO", "PZA", "EMB", "STR", "RIF", "INH"], color = 'black', fc = 'none', header = 'TB-drug data (HBW)')
plt.savefig(f"tb_drug_plot_hyperHMM.png", dpi=900)
plt.show()

bubble_plot(10, txt2matrix(f"Plot files/hypertraps_tb.txt"), names = ["AMI", "CAP", "MOX", "OFL", "PRO", "PZA", "EMB", "STR", "RIF", "INH"], color = 'blue', header = f'TB-drug data (HyperTraPS)')
plt.grid()
plt.savefig(f"tb_drug_plot_hypertraps.png", dpi=900)
plt.show()




mean = ["mean_simple2_L5.txt", "mean_simple2_L7.txt", "mean_simple2_L9.txt", "mean_double2_L5.txt" , "mean_double2_L7.txt", "mean_double2_L9.txt"]
sd = ["sd_simple2_L5.txt", "sd_simple2_L9.txt", "sd_simple2_L9.txt", "sd_double2_L5.txt", "sd_double2_L9.txt", "sd_double2_L9.txt"]

for i in range(len(mean)):
    if i < 3:
        bubble_plot(5+i*2, txt2matrix("Plot files/" + mean[i]), names = list(range(1, 5+i*2+1)), color = 'dodgerblue', header = f'Simple case samples (HBW)')
        bubble_plot(5+i*2, txt2matrix("Plot files/" + sd[i]), names = list(range(1, 5+i*2+1)), color = 'black', fc = 'none', header = f'Simple case samples (HBW)')
        plt.savefig(f"single_L{5+i*2}_hyperHMM.png", dpi=900)
        plt.show()
    else:
        bubble_plot(5+(i-3)*2, txt2matrix("Plot files/" + mean[i]), names = list(range(1, 5+(i-3)*2+1)), color = 'dodgerblue', header = f'Double case samples (HBW)')
        bubble_plot(5+(i-3)*2, txt2matrix("Plot files/" + sd[i]), names = list(range(1, 5+(i-3)*2+1)), color = 'black', fc = 'none', header = f'Double case samples (HBW)')
        plt.savefig(f"double_L{5+(i-3)*2}_hyperHMM.png", dpi=900)
        plt.show()




def hypertraps_results2matrix(csv_name, L):
    data = pd.read_csv(csv_name)
    matrix = np.empty((L,L))
    for i in range(L):
        for j in range(L):
            matrix[i,j] = data["Prob"].to_numpy()[i*L+j]
        
    return matrix


hyp_L = [5,7,9]
for i in hyp_L:
    bubble_plot(i, hypertraps_results2matrix(f"Plot files/cross-{i}.csv", i),names = list(range(1, i+1)), header= "Double case samples (HyperTraPS)")
    plt.grid()
    plt.savefig(f"double_L{i}_hypertraps.png", dpi=900)
    plt.show()




