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
    bubble_plot(5, txt2matrix(f"mean_simple{i}.txt"), names = [1,2,3,4,5], color = 'dodgerblue', header = f'Simple case {i*4} samples (HBW)')
    bubble_plot(5, txt2matrix(f"sd_simple{i}.txt"), names = [1,2,3,4,5], color = 'black', fc = 'none', header = f'Simple case {i*4} samples (HBW)')
    plt.show()
    bubble_plot(5, txt2matrix(f"hypertraps_single_L5_{i}.txt"), color = 'blue', header = f'HyperTraPS {i*4} samples')
    plt.grid()
    plt.show()





#Ovarian cancer plotting:

bubble_plot(7, txt2matrix("mean_ovarian.txt"), names = ['Xp-', '1q+', '8p-', '4q-', '5q-', '3q+', '8q+'], color = 'dodgerblue', header = 'Ovarian cancer data (HBW)')
bubble_plot(7, txt2matrix("sd_ovarian.txt"), names = ['Xp-', '1q+', '8p-', '4q-', '5q-', '3q+', '8q+'], color = 'black', fc = 'none', header = 'Ovarian cancer data (HBW)')
plt.show()

bubble_plot(7, txt2matrix(f"hypertraps_ovarian.txt"),  names = ['Xp-', '1q+', '8p-', '4q-', '5q-', '3q+', '8q+'], color = 'blue', header = f'Ovarian cancer data (HyperTraPS)')
plt.grid()
plt.show()




#TB-drug plotting:

bubble_plot(10, txt2matrix("mean_tb_drug.txt"), names = ["AMI", "CAP", "MOX", "OFL", "PRO", "PZA", "EMB", "STR", "RIF", "INH"], color = 'dodgerblue', header = 'TB-drub data (HBW)')
bubble_plot(10, txt2matrix("sd_tb_drug.txt"), names = ["AMI", "CAP", "MOX", "OFL", "PRO", "PZA", "EMB", "STR", "RIF", "INH"], color = 'black', fc = 'none', header = 'TB-drug data (HBW)')
plt.show()

bubble_plot(10, txt2matrix(f"hypertraps_tb.txt"), names = ["AMI", "CAP", "MOX", "OFL", "PRO", "PZA", "EMB", "STR", "RIF", "INH"], color = 'blue', header = f'TB-drug data (HyperTraPS)')
plt.grid()
plt.show()




#Timing test plotting:
mean = ["mean_simple2_L5.txt", "mean_simple2_L7.txt", "mean_simple2_L9.txt", "mean_double2_L5.txt" , "mean_double2_L7.txt", "mean_double2_L9.txt"]
sd = ["sd_simple2_L5.txt", "sd_simple2_L9.txt", "sd_simple2_L9.txt", "sd_double2_L5.txt", "sd_double2_L9.txt", "sd_double2_L9.txt"]

for i in range(len(mean)):
    if i < 3:
        bubble_plot(5+i*2, txt2matrix(mean[i]), names = list(range(1, 5+i*2+1)), color = 'dodgerblue', header = f'Simple case samples (HBW)')
        bubble_plot(5+i*2, txt2matrix(sd[i]), names = list(range(1, 5+i*2+1)), color = 'black', fc = 'none', header = f'Simple case samples (HBW)')
        plt.show()
    else:
        bubble_plot(5+(i-3)*2, txt2matrix(mean[i]), names = list(range(1, 5+(i-3)*2+1)), color = 'dodgerblue', header = f'Double case samples (HBW)')
        bubble_plot(5+(i-3)*2, txt2matrix(sd[i]), names = list(range(1, 5+(i-3)*2+1)), color = 'black', fc = 'none', header = f'Double case samples (HBW)')
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
    bubble_plot(i, hypertraps_results2matrix(f"cross-{i}.csv", i),names = list(range(1, i+1)), header= "Double case samples (HyperTraPS)")
    plt.grid()
    plt.show()








#Timing results

time_hypertraps_01 = [168, 950, 1308, 305, 1081, 2521]
time_hypertraps_001 = [961, 3376, 8896, 1181, 5914, 16978]
time_abw = [0.140, 1.59, 33.6, 0.397, 3.18, 67.5]


x_lab = ["S L=5", "S L=7", "S L=9", "D L=5", "D L=7", "D L=9"]

N = 6
ind = np.arange(N) 
width = 0.25
  
bar1 = plt.bar(ind, np.log(time_abw), width, color = 'r')
  
bar2 = plt.bar(ind+width, np.log(time_hypertraps_001), width, color='g')
  
bar3 = plt.bar(ind+width*2, np.log(time_hypertraps_01), width, color = 'b')
  
plt.xlabel("Test", fontsize=18)
plt.ylabel('log(time) [s]', fontsize=18)
plt.title("Timing results", fontsize=20)
  
plt.xticks(ind+width,x_lab)
plt.legend( (bar1, bar2, bar3), ('HBW', 'HyperTraPS 0.001', 'HyperTraPS 0.01') )
plt.show()


#Timing results for simple case and higher L's

x = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
y = [0.0018, 0.0036, 0.0073, 0.022, 0.035, 0.089, 0.24, 0.72, 1.88, 5.54, 13.46, 33.74, 85.58, 236.14, 577.29, 1447.79]


plt.plot(x,y, "-o")
  
plt.xlabel("L", fontsize=18)
plt.ylabel('time [s]', fontsize=18)
plt.title("Timing results", fontsize=20)  
plt.show()