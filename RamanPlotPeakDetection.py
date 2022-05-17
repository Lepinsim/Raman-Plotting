### Plots spectra in every directories at 'filepath' string address
### Does baseline substraction, normalization, peakdetection and plotting

import scipy as scipy
from scipy import optimize
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from sklearn.preprocessing import StandardScaler, MinMaxScaler, Normalizer
from sklearn.preprocessing import RobustScaler

# from scipy.interpolate import splev

import math
import os
import seaborn as sns
# from matplotlib.colors import ListedColormap
# from scipy.signal import find_peaks

from scipy import sparse
from scipy.sparse.linalg import spsolve
# from scipy import signal
# from adjustText import adjust_text
# from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks



def baseline_als(y, lam=100000, p=0.1, niter=10):
  L = len(y)
  D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
  w = np.ones(L)
  for i in range(niter):
    W = sparse.spdiags(w, 0, L, L)
    Z = W + lam * D.dot(D.transpose())
    z = spsolve(Z, w*y)
    w = p * (y > z) + (1-p) * (y < z)
  return z

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w



#####################################################################################################################""
#####################################################################################################################""
#####################################################################################################################""
#####################################################################################################################""

 ## Initialization
# Working folder

path = r'C:\Users\Simon\surfdrive\Experiment\Mix FeSO4\feso4mix_Raman\1analysis\New folder'

Figscale =0.2 
nbPeakshow = 15
MaxHeight = 5
SpectraSpacing = 15 

#####################################################################################################################""
#####################################################################################################################""
#####################################################################################################################""

    # Start looping over folders containing spectra

directories = []
print(os.walk(path))
for subdir, dirs, files in os.walk(path):
    # directories = []
    print(subdir,'end')

    names = []
    for file in files:
        # print(os.path.join(subdir, file))
        filepath = os.path.join(subdir, file)
        if filepath.endswith('.csv'):
          names.append(filepath)
    directories.append(names)
    
    
i=1
nSpectra = len(names)

# inialization empty containers
dfPeaks = pd.DataFrame()
dfPeaksPrint = pd.DataFrame()
emptyCol = np.linspace(0,14,15)



sns.set_style("white")
sns.set_style("ticks")


snscolor=  sns.diverging_palette(200, 15, s=75, l=60,n=nSpectra + 1, center="dark")
# snscolor=  sns.color_palette("tab10")
# snscolor = ['#0000ffff','#7d7dffff','#ff7d7dff','#ff0000ff']


#####################################################################################################################""
#####################################################################################################################""
#####################################################################################################################""
 ### # Starting loop loading and plotting each curve in each directory 
#Iterating over folders
for names in directories:

    # creating an empty plot
  fig, ax1 = plt.subplots()
    #Iterating over all files
  for j in range (len(names)):
    #Plotting only one curve
  # for j in range(1):

        # Initializing variables
        name = names[j]
        # print(name)
        dfa = pd.DataFrame()
         # Interval between lines ----------------------------------------------------------
        deltaY = SpectraSpacing*(-0.3-1/nSpectra)*i

        #Number from 0 to 1 scaling figure frame and fontsize 


        ##Loading data/
        dfb = pd.read_csv(name,header=None, skiprows=19, encoding='cp1252')


         ### Preprocessing
         # Baseline removal
         #Loads x axis
        dfa.loc[:,0] = dfb.loc[:,0]
         #Loads Y axis and substract a fitted baseline in dfa
        dfa.loc[:,1] = dfb.iloc[:,1]-baseline_als(dfb.iloc[:,1],100000,0.001)

    
         # Scale intensities to [0,1]for normalization
        dfanorm = StandardScaler().fit_transform(dfa.values)[:,1]

        #actual plotting function plotting x and normalized Y
        p, = ax1.plot(dfa.iloc[:,0],dfanorm-deltaY,
            linewidth=1,
            color=snscolor[i],
            label=name[82:-4])

         
  #
        ###   Detection of the peaks

        peaks2, p2H_dict = find_peaks(dfanorm, distance=10, prominence=0.5, height=0.15, width=1) #--------------------------------------

        #data managing (not interesting)
        p2H_list = list(p2H_dict.items())
        p2H_array = np.array(p2H_list[0][1])
        peaks2b = np.array(list(zip(peaks2,p2H_array)))
        
        # Peak list sort in ascending order
        peaks3 = peaks2b[peaks2b[:, 1].argsort()]
        print(peaks3)

      ###########################k## PLOT DESIGN #################
          #### write each peak

        # possibility to save peaklist        
       # np.savetxt(r'C:\Users\Anne\Desktop\Stage uva\april\22-04\Raman\avec EtOH\fichier spectre\25Fe\1 ' + str(i) + '.csv', peaks3, delimiter=",")
        

        
        for peak in peaks3[-nbPeakshow:,0]:
          
          PeakLblY = dfanorm[int(peak)]
          if PeakLblY > MaxHeight:
              PeakLblY = MaxHeight
          # ax1.annotate(int(dfa.iloc[int(peak),0]),(dfa.iloc[int(peak),0], PeakLblY-deltaY))
          ax1.annotate(int(dfa.iloc[int(peak),0]),
                       (dfa.iloc[int(peak),0], dfanorm[int(peak)]-deltaY),
                       color=snscolor[i], 
                       fontsize=40*Figscale)
       
        # plt.savefig(r'C:\Users\Simon\surfdrive\Experiment\Mix FeSO4\FeSO4_NaCl early2022\Mar22')

        i +=1

  print(nSpectra)
  fontsizeAxe = 40*Figscale
  fig.set_figheight(40*Figscale)
  fig.set_figwidth(30*Figscale)  
  plt.xticks(fontsize=fontsizeAxe)
  plt.yticks(fontsize=fontsizeAxe)
  # plt.xticks()
  # plt.yticks()
  #plt.xlim(None  ,1500)
  # plt.ylim(None,6)


  handles, labels = ax1.get_legend_handles_labels()
  ax1.legend(reversed(handles), reversed(labels), loc='upper right', fontsize=fontsizeAxe*0.5)


  plt.xlabel('Raman Shift (rel. cm$^{-1}$)', fontsize=fontsizeAxe+2, fontweight='bold')
  plt.ylabel('Intensity (a. u.)', fontsize=fontsizeAxe+2, fontweight='bold')

plt.show()




