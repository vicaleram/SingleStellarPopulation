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
    idx1 = i.index('_03_')
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
    #omit zero values in the dataframe and other weird stars
    data.drop(data.query("`log_Teff`<4.3 and `log_L`>4.0").index)
    data = data.loc[data['log_Teff']>3.90]
    datafiles[str(ages[j])] = data
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

#do the procedure of determining the apparent magnitude of each star to find their minimum detectable B field
##Calculate bolumetric corrections
##define the coeffecients from the table of Torres(2010)
c0 = [-0.190537291496456e5, 0.155144866764412e5, 0.421278819301717e4, 0.381476328422343e3]
c1 = [-0.370510203809015e5, 0.385672629965804e5, -0.150651486316025e5, 0.261724637119416e4, -0.170623810323864e3]
c2 = [ -0.118115450538963e6, 0.137145973583929e6, -0.636233812100225e5, 0.147412923562646e5, -0.170587278406872e4, 0.788731721804990e2]
#define the correction for the sun and its apparent magntiude
v_sun = -26.76
bc_sun = -0.06

#DEFINE A DISTANCE FOR THE STAR CLUSTER IN PARSECS, THIS IS ARBITRARY
distance = 400   #pc

for age in ages:
    cluster = datafiles[str(age)]
    cluster['BC'] = 1.0
    logT_0 = cluster.loc[cluster['log_Teff']<=3.70,'log_Teff']
    logT_1 = cluster[cluster['log_Teff'].between(3.70,3.90)]
    #print(cluster[cluster['log_Teff'].between(3.70,3.90)])
    logT_2 = cluster.loc[cluster['log_Teff']>=3.90,'log_Teff']
    bc0 = c0[0] +c0[1]*logT_0 +c0[2]*(logT_0)**2 + c0[3]*(logT_0)**3
    bc1 = c1[0] +c1[1]*logT_1 +c1[2]*(logT_1)**2 + c1[3]*(logT_1)**3 +c1[4]*(logT_1)**4
    bc2 = c2[0] +c2[1]*logT_2 +c2[2]*(logT_2)**2 + c2[3]*(logT_2)**3 +c2[4]*(logT_2)**4 + c2[5]*(logT_2)**5
    cluster.loc[cluster['log_Teff']<3.70, 'BC'] = bc0
    cluster.loc[cluster['log_Teff'].between(3.70,3.90), 'BC'] = bc1
    cluster.loc[cluster['log_Teff']>=3.90,'BC'] = bc2
    #check that both the cluster DF and the Dic are the same
    #print(cluster['BC'].iloc[0])
    #print(datafiles[str(age)]['BC'].iloc[0], 'Dic')
    ### calcultae the Absolute Visual magntiude for all stars 
    M_v = -2.5*cluster['log_L'] + v_sun +31.572-(cluster['BC']- bc_sun)
    cluster['M_V'] = M_v
    
    ###calculate apparent visul magnitude given the distance defined above
    m_v = cluster['M_V'] + 5*np.log10(distance/10)
    cluster['m_V'] = m_v
    cluster.drop(cluster[cluster['m_V']>20].index)

'''Assign inclinations to all the stars and keep them constant'''
#Assign an inclination for the initial data file so that stars keep constant inclination
ini_data= datafiles[str(ages[0])]
c = np.random.rand(len(ini_data))
i_val = np.arccos(1 - c)
ini_data['i']=i_val
ini_data['vsini']= ini_data['surf_avg_v_rot']*np.sin(ini_data['i'])

for k in range(1,len(ages)):
    cluster = datafiles[str(ages[k])].copy()
    #assign their respective inclination and calculate vsini
    cluster = pd.merge(cluster, ini_data[['#Star_ID','i']], on='#Star_ID').copy()
    cluster['vsini'] = cluster['surf_avg_v_rot']*np.sin(cluster['i'])
    datafiles[str(ages[k])] = cluster

'''This part is used to compute the minimum B field that can be detected using MOBSTER values
'''
#Define the exposure time. IT IS ARBITRARY, CHOOSE WHICEVER VALUE IS REASONABLE
#define exposure time in seconds
t_exp = 7200 # 2hours in seconds
s_exp = t_exp /4  #sub exposure time, simplye T_exp/4
#We need a linear fit to obtain the LSD gain as a function of temperature. Full explanation is in original notebook
slope = -12.388115933619902
intercept = 61.67535680103128
mag_inc = []
for age in ages:
    cluster = datafiles[str(age)].copy()
    snr_bin = 430*10**((8.4-cluster['m_V'])/5)*np.sqrt(4*s_exp/3200) #Signal To Noise Ratio bin 
    #print(snr_bin)
    cluster['SNR_b'] = snr_bin
    snr_18_a = cluster['SNR_b'][cluster['vsini']<=18]
    snr_18_b = cluster['SNR_b'][cluster['vsini']>18]*np.sqrt(cluster['vsini'][cluster['vsini']>18]/18)
    cluster['SNR_1.8'] = np.nan
    cluster.loc[cluster['vsini']<=18, 'SNR_1.8'] = snr_18_a
    cluster.loc[cluster['vsini']>18, 'SNR_1.8'] = snr_18_b
    #calculate the SNR for the LSD gain 
    lsd_gain = intercept + slope*cluster['log_Teff']
    snr = cluster['SNR_1.8']*lsd_gain
    cluster['SNR']= snr
    #print(snr)
    #estimate the detectable magnetic field
    b_min_a = (200+ 289*cluster['vsini'][cluster['vsini']<= 40])/cluster['SNR'][cluster['vsini']<= 40]
    #print(b_min_a)
    b_min_b = (-32243 + 1077.9*cluster['vsini'][cluster['vsini']> 40])/cluster['SNR'][cluster['vsini']>40]
    #print(b_min_b)
    cluster['B_min'] = np.nan
    cluster.loc[cluster['vsini']<=40, 'B_min'] = b_min_a*100
    cluster.loc[cluster['vsini']>40, 'B_min'] = b_min_b*100
    nondete_n = len(cluster[cluster['B_min']>cluster['Beq']])
    #add star that start with 0 b field
    nonmag_n = len(cluster[cluster['Binit']<= 0.1])
    mag_per = (nondete_n+nonmag_n)/len(cluster)
    mag_inc.append(mag_per)
    datafiles[str(age)] = cluster


#Make animation of evolution of minimum magentic field of the cluster
#make animation of evolution of Stellar rotation with age
fig, ax = plt.subplots(1, 1, figsize=(10,8))
ini_data = datafiles[str(ages[0])]
##plot the initial distribution of Bmin
#print(ini_data['B_min'])
ax.hist(ini_data['B_min'],30, color= 'b', edgecolor = 'k');
ax.set_xlabel(r'$B_\mathrm{min}$ (G)')
ax.set_ylabel('N')
ax.text(np.max(ini_data['B_min'])-100, 250, r'Age $={:.2e}$yr'.format(ages[0]), size=20, bbox =dict(boxstyle='round', facecolor = 'white'))
ax.set_xlim(0,np.max(ini_data['B_min']))
ax.set_ylim(0,850)
ax.text(1250, 250, r'Age $={:.2e}$yr'.format(ages[0]), size=20, bbox =dict(boxstyle='round', facecolor = 'white'))

def animate_Bmin(i):
    ax.cla()
    current = datafiles[str(ages[i])]
    ax.hist(current['B_min'],30, color= 'b', edgecolor = 'k');
    ax.set_xlabel(r'$B_\mathrm{min}$ (G)')
    ax.set_ylabel('N')
    ax.set_xlim(0,np.max(ini_data['B_min'])-100)
    ax.set_ylim(0,850)
    ax.text(np.max(ini_data['B_min'])-100, 250, r'Age $={:.2e}$yr'.format(ages[i]), size=20, bbox =dict(boxstyle='round', facecolor = 'white'))
    #ax.set_xlim(0,1000)
    #ax.set_ylim(0,700)
ani = animation.FuncAnimation(
    fig, animate_Bmin, interval=len(ages), save_count=len(ages), frames = len(ages))
plt.draw()
plt.show()

## To save the animation
writergif = animation.PillowWriter(fps=1)
ani.save(figpath+'Bmin_'+new_name+'.gif', writer=writergif)



#make HR of detectable stars
fig , ax = plt.subplots(figsize=(10,10))
ax.scatter(ini_data.loc[ini_data['Beq']< ini_data['B_min'], 'log_Teff'],ini_data.loc[ini_data['Beq']< ini_data['B_min'], 'log_L'], color = 'r', s= 0.2*(10**ini_data.loc[ini_data['Beq']< ini_data['B_min'],'log_R'])**2, label = 'Non-detectable', alpha = 0.7)
ax.scatter(ini_data.loc[ini_data['Beq']>ini_data['B_min'], 'log_Teff'], ini_data.loc[ini_data['Beq']>ini_data['B_min'], 'log_L'],color = 'g', s= 0.2*(10**ini_data.loc[ini_data['Beq']> ini_data['B_min'],'log_R'])**2, label = 'Detectable', alpha = 0.7)
ax.scatter(ini_data.loc[ini_data['Binit']<= 0.1, 'log_Teff'],ini_data.loc[ini_data['Binit']<= 0.1, 'log_L'], color = 'b', s= 0.2*(10**ini_data.loc[ini_data['Binit']<= 0.1,'log_R'])**2, label = 'Non-Magnetic', alpha = 0.7)
ax.set_ylim(1.8,6.1)
ax.set_xlim(4.05,4.68)
ax.set_xlabel(r'$\log_{10}{T_{eff}}$')
ax.set_ylabel(r'$\log_{10}{L}$')
ax.invert_xaxis()
#ax.set_title('HR diagram for {:.1e} stars cluster'.format(n_stars))
ax.text(4.3, 5.0, r'Age $={:.2e}$yr'.format(ages[0]),weight = 'bold', size=20, bbox =dict(boxstyle='round', facecolor = 'white'))
ax.text(4.55, 2.5, r'$\%$ Incidence'.format(mag_inc[0]), size=20, bbox =dict(boxstyle='round', facecolor = 'white'))
ax.legend(fontsize = 20)
#plot MESA models
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

def animate_HRdet(i):
    ax.cla()
    current = datafiles[str(ages[i])]
    ax.scatter(current.loc[current['Beq']< current['B_min'], 'log_Teff'],current.loc[current['Beq']< current['B_min'], 'log_L'], color = 'r', s= 0.2*(10**current.loc[current['Beq']< current['B_min'],'log_R'])**2, label = 'Non-detectable', alpha = 0.7)
    ax.scatter(current.loc[current['Beq']>current['B_min'], 'log_Teff'], current.loc[current['Beq']>current['B_min'], 'log_L'],color = 'g', s= 0.2*(10**current.loc[current['Beq']> current['B_min'],'log_R'])**2, label = 'Detectable', alpha = 0.7)
    ax.scatter(current.loc[current['Binit']<= 0.1, 'log_Teff'],current.loc[current['Binit']<= 0.1, 'log_L'], color = 'b', s= 0.2*(10**current.loc[current['Binit']<= 0.1,'log_R'])**2, label = 'Non-Magnetic', alpha = 0.7)
    ax.set_ylim(1.8,6.1)
    ax.set_xlim(4.05,4.68)
    ax.invert_xaxis()
    ax.set_xlabel(r'$\log_{10}{T_{eff}}$')
    ax.set_ylabel(r'$\log_{10}{L}$')
    #ax.set_title('HR diagram for {:.1e} stars cluster'.format(n_stars))
    ax.text(4.3, 5.0, r'Age $={:.2e}$yr'.format(ages[i]),weight='bold', size=20, bbox =dict(boxstyle='round', facecolor = 'white'))
    ax.text(4.55, 2.5, r'${:.2f}\%$ Incidence'.format((1-mag_inc[i])*100), size=20, bbox =dict(boxstyle='round', facecolor = 'white'))
    ax.legend(fontsize = 20)
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

ani = animation.FuncAnimation(
    fig, animate_HRdet, interval=len(ages), save_count=len(ages), frames = len(ages))
plt.draw()
plt.show()

#print(mag_inc)
## To save the animation
writergif = animation.PillowWriter(fps=8)
ani.save(figpath+'detectability_'+new_name+'.gif', writer=writergif)

#plot the detectabiliyt as a function of time
fig , ax = plt.subplots(figsize=(8,8))
ax.plot(ages,(1-np.array(mag_inc))*100, color = 'b')
ax.set_xlabel('Cluster age [yr]')
ax.set_ylabel('Obsevational incidence [%]')
ax.set_title('Time evolution magnetic incidence')
plt.tight_layout()
plt.savefig(figpath+'incidence'+new_name+'.pdf', format = 'pdf')

