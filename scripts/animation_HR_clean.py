# imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl
#plt.rcParams['animation.ffmpeg_path'] = '/home/viclrd/.local/lib/python3.9/site-packages/ffmpeg'
from matplotlib import cm
import sys
import os
import re
#print(os.path.expanduser('~'))

#make better plots
plt.style.use(['science','notebook','grid'])
plt.rcParams.update({"font.size" : 26})
plt.rcParams.update({"axes.labelsize" : 26})
plt.rcParams.update({"xtick.labelsize" : 22})
plt.rcParams.update({"ytick.labelsize" : 22})
plt.rcParams.update({"axes.titlesize" : 26})
plt.rcParams.update({"figure.figsize": (10,8)})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'


#specify the path
path = sys.argv[1]
name = sys.argv[1]
#separate the string in case of giving full path
separate = name.split('/')
new_name = separate[-2]
#path of images
figpath = os.path.expanduser('~')+'/figures/'
#separate the string in case of giving full path
separate = path.split('/')
folder = separate[-2]

#open all the data files inside the folder passed and store them
names = [filename for filename in os.listdir(path)]

#find the ages of each file
ages = []
for i in names:
    idx1 = i.index('_01_')
    idx2 = i.index('yrs')
    temp = float(i[idx1 +len('_00_'):idx2])
    ages.append(temp)

#read each data file and make a dataframe
#create dictionary to store them
datafiles = {}
j=0
for filename in names:
    data = pd.read_csv(path+filename, sep='\t', header = 16)
    data.columns = data.columns.str.replace(' ','')
    #omit zero values in the dataframe
    #data.drop(data.query("`log_Teff`<4.3 and `log_L`>4.0").index)
    #data = data.loc[data['log_Teff']>3.90]
    datafiles[str(ages[j])] = data
    del(data)
    j = j +1

ages =np.sort( np.array(ages))


# read the MESA models that have no magnetic field to overplot
grid = 'NOMAGZ14Mix1/0/'
mesa_models = np.concatenate((np.arange(3.0,26.0,2),[30,40,50,60]))
mesa_names = []
for i in mesa_models:
    #print( '{:.3}'.format(i))
    temp = '{:.3}'.format(i)
    mesa_names.append(temp)

#make animation of HR diagram at each age
fig, ax = plt.subplots(figsize=(11,11))
#first make it the initial age
lower = np.linspace(0,1000,5)
middle = np.linspace(1500,10000,18)
upper = np.array([15e3,2e4,3e4,5e4])
bounds = np.concatenate((lower,middle,upper),axis = None)
N = len(bounds)
cmap = plt.get_cmap('plasma', N)
ini_data = datafiles[str(ages[0])]
n_stars = ini_data.shape[0]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

ax.set_ylim(1.8,6.1)
ax.set_xlim(4.05,4.68)
ax.set_xlabel(r'$\log_{10}{T_{eff}}$')
ax.set_ylabel(r'$\log_{10}{L}$')
ax.set_title('HR diagram for {:.1e} stars cluster'.format(n_stars))
ax.invert_xaxis()
for j in mesa_names:
    mesa = pd.DataFrame(np.genfromtxt(grid+j+'/LOGS/out.data', skip_header=5, names=True))
    #get the zams
    h0 = np.max(mesa['center_h1'])
    h_zams = h0 - h0*0.3/100
    mesa = mesa[mesa['center_h1'] <= h_zams]
    ax.plot(mesa['log_Teff'],mesa['log_L'], ls = '--', color = 'k')
    if float(j) < 30:
        ax.text(np.max(mesa['log_Teff'])+0.07,np.min(mesa['log_L'])-0.01, r'$'+j+'M_\odot$',size=16)
    else:
        ax.text(np.min(mesa['log_Teff'])-0.02,np.max(mesa['log_L'])-0.08, r'$'+j+'M_\odot$',size=16)
ax.text(4.6, 2.5, r'Age $={:.2e}$yr'.format(ages[0]), size=20, bbox =dict(boxstyle='round', facecolor = 'white'))
area =0.2* (10**ini_data['log_R'])**2
ban = ax.scatter(ini_data['log_Teff'],ini_data['log_L'],\
                 c =ini_data['Beq'], cmap=cmap,s=area, norm=norm,alpha = 0.7)

fig.colorbar(ban, ticks = bounds,label=r'$B (\mathrm{G})$')

#define the animation function
def animate(i):
    plt.cla()
    current = datafiles[str(ages[i])]
    n_stars = current.shape[0]
    ax.set_ylim(1.8,6.1)
    ax.set_xlim(4.05,4.68)
    ax.set_xlabel(r'$\log_{10}{T_{eff}}$')
    ax.set_ylabel(r'$\log_{10}{L}$')
    ax.set_title('HR diagram for {:.1e} stars cluster'.format(n_stars))
    ax.invert_xaxis()
    for j in mesa_names:
        mesa = pd.DataFrame(np.genfromtxt(grid+j+'/LOGS/out.data', skip_header=5, names=True))
        h0 = np.max(mesa['center_h1'])
        h_zams = h0 - h0*0.3/100
        mesa = mesa[mesa['center_h1'] <= h_zams]
        ax.plot(mesa['log_Teff'],mesa['log_L'], ls = '--', color = 'k')
        if float(j) < 30:
            ax.text(np.max(mesa['log_Teff'])+0.07,np.min(mesa['log_L'])-0.01, r'$'+j+'M_\odot$',size=16)
        else:
            ax.text(np.min(mesa['log_Teff'])-0.02,np.max(mesa['log_L'])-0.08, r'$'+j+'M_\odot$',size=16)
    ax.text(4.6, 2.5, r'Age $={:.2e}$yr'.format(ages[i]), size=20, bbox =dict(boxstyle='round', facecolor = 'white'))
    area =0.2* (10**current['log_R'])**2
    ax.scatter(current['log_Teff'],current['log_L'],\
                 c =current['Beq'], cmap=cmap,s=area, norm=norm,alpha = 0.7)
    plt.tight_layout()

ani = animation.FuncAnimation(
    fig, animate, interval=len(ages), save_count=len(ages), frames = len(ages))
plt.draw()
plt.show()

# To save the animation, use e.g.
#
writergif = animation.PillowWriter(fps=8)
ani.save(figpath+'hr_'+new_name+'.gif', writer=writergif)

#Do the same but now for the surface rotation of the stars
fig, ax = plt.subplots(figsize=(11,11))
n_stars = ini_data.shape[0]
cmap = plt.get_cmap('viridis')
ax.set_ylim(1.8,6.1)
ax.set_xlim(4.05,4.68)
ax.set_xlabel(r'$\log_{10}{T_{eff}}$')
ax.set_ylabel(r'$\log_{10}{L}$')
ax.set_title('HR diagram for {:.1e} stars cluster'.format(n_stars))
ax.invert_xaxis()
for j in mesa_names:
    mesa = pd.DataFrame(np.genfromtxt(grid+j+'/LOGS/out.data', skip_header=5, names=True))
    h0 = np.max(mesa['center_h1'])
    h_zams = h0 - h0*0.3/100
    mesa = mesa[mesa['center_h1'] <= h_zams]
    ax.plot(mesa['log_Teff'],mesa['log_L'], ls = '--', color = 'k')
    if float(j) < 30:
        ax.text(np.max(mesa['log_Teff'])+0.06,np.min(mesa['log_L'])-0.005, r'$'+j+'M_\odot$',size=16)
    else:
        ax.text(np.min(mesa['log_Teff'])-0.02,np.max(mesa['log_L'])-0.08, r'$'+j+'M_\odot$',size=16)
ax.text(4.6, 2.5, r'Age $={:.2e}$yr'.format(ages[0]), size=20, bbox =dict(boxstyle='round', facecolor = 'white'))
area =0.2* (10**ini_data['log_R'])**2
ban = ax.scatter(ini_data['log_Teff'],ini_data['log_L'],\
                 c =ini_data['surf_avg_v_rot'], cmap=cmap,s=area,alpha = 0.7)

fig.colorbar(ban,label=r'Surface $v_\mathrm{avg}$ (km/s)')

#define the animation function
def animate_v(i):
    plt.cla()
    current = datafiles[str(ages[i])]
    n_stars = current.shape[0]
    ax.set_ylim(1.8,6.1)
    ax.set_xlim(4.05,4.68)
    ax.set_xlabel(r'$\log_{10}{T_{eff}}$')
    ax.set_ylabel(r'$\log_{10}{L}$')
    ax.set_title('HR diagram for {:.1e} stars cluster'.format(n_stars))
    ax.invert_xaxis()
    for j in mesa_names:
        mesa = pd.DataFrame(np.genfromtxt(grid+j+'/LOGS/out.data', skip_header=5, names=True))
        h0 = np.max(mesa['center_h1'])
        h_zams = h0 - h0*0.3/100
        mesa = mesa[mesa['center_h1'] <= h_zams]
        ax.plot(mesa['log_Teff'],mesa['log_L'], ls = '--', color = 'k')
        if float(j) < 30:
            ax.text(np.max(mesa['log_Teff'])+0.06,np.min(mesa['log_L'])-0.005, r'$'+j+'M_\odot$',size=16)
        else:
            ax.text(np.min(mesa['log_Teff'])-0.02,np.max(mesa['log_L'])-0.08, r'$'+j+'M_\odot$',size=16)
    ax.text(4.6, 2.5, r'Age $={:.2e}$yr'.format(ages[i]), size=20, bbox =dict(boxstyle='round', facecolor = 'white'))
    area =0.2* (10**current['log_R'])**2
    ax.scatter(current['log_Teff'],current['log_L'],\
                 c =current['surf_avg_v_rot'], cmap=cmap,s=area,alpha = 0.7)
    plt.tight_layout()

ani = animation.FuncAnimation(
    fig, animate_v, interval=len(ages), save_count=len(ages), frames = len(ages))
plt.draw()
plt.show()

# To save the animation, use e.g.
#
writergif = animation.PillowWriter(fps=8)
ani.save(figpath+'hr_v_'+new_name+'.gif', writer=writergif)
