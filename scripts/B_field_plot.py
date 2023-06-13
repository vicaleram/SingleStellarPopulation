#make a script to plot the stars from a simulation in an HR diagram

#imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

#make prettier plots
plt.style.use(['science','notebook','grid'])
plt.rcParams.update({"font.size" : 26})
plt.rcParams.update({"axes.labelsize" : 24})
plt.rcParams.update({"xtick.labelsize" : 22})
plt.rcParams.update({"ytick.labelsize" : 22})
plt.rcParams.update({"axes.titlesize" : 26})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

#specify the path
path = './'
name = sys.argv[1]
figpath = '../figures/'

#read the data file and make a dataframe
data = pd.read_csv(name, sep='\t', header = 16)
data.columns = data.columns.str.replace(' ', '')

#make the B field distribution
import matplotlib as mpl
from matplotlib import cm
fig, ax = plt.subplots(figsize=(10,10))

n_stars = data.shape[0]

#ax.set_ylim(2.6,np.nanmax(data['log_L'])+0.1)
#ax.set_ylim(2.55,np.sort(data['log_L'])[-2]+0.1)
#ax.set_xlim(4.235,np.nanmax(data['log_Teff'])+0.05)
ax.set_title('B field distribution for {:.1e} stars'.format(n_stars))
ax.hist(data['Binit'], 50, color = 'b', edgecolor = 'k')
ax.set_xlabel(r'$B(G)$');
ax.set_ylabel(r'Number of Stars')
plt.tight_layout()
new_name = name[0:-4]
plt.savefig(figpath+'b_f'+new_name+'.pdf', format = 'pdf')

