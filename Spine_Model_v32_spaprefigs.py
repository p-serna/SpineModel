#!/usr/bin/env python
# coding: utf-8

# # Figures for the paper

# In[1]:


from scipy.integrate import quad
from numpy import *
from matplotlib.pylab import *
import matplotlib as mpl
import os
import sys

# changing parameters

try:
    assert len(sys.argv)==2, 'No correct number of arguments'
    from data_Fullset import SCxL23 as exppar
    print('No shrinkage correction')
    shflag = ''
except:
    assert len(sys.argv)==3, 'No correct number of arguments'
    shrinkage = float(sys.argv[2])
    if shrinkage > 0.0:
        print('We use corrected shrinkage')
        shflag = 'sh'    
    else:
        from data_Fullset import SCxL23 as exppar
        print('No shrinkage correction')
        shflag = ''    

condition = sys.argv[1]
#condition = '1000_80_65'

folderstore = './fconditions/spatial/'

    

# In[2]:


import PIL.Image as Image
import pickle


# In[3]:


def plot_trace(data,t,ax=None,c='C0',band= None):
    if ax is None:
        ax = gca()
    vtracso = data*1.0
    vtracso.sort(axis=1)
    ax.plot(t,vtracso[:,250],c)
    if band == 0:
        pass
    elif band == 1:
        ax.fill_between(t,vtracso[:,79],vtracso[:,421],color=c,alpha=0.5)
    elif band ==2:
        ax.fill_between(t,vtracso[:,5],vtracso[:,495],color=c,alpha=0.2)
    else:
        ax.fill_between(t,vtracso[:,5],vtracso[:,495],color=c,alpha=0.2)
        ax.fill_between(t,vtracso[:,79],vtracso[:,421],color=c,alpha=0.5)
    return(ax)


# In[4]:


from matplotlib.patches import ConnectionPatch


# In[5]:


def getint(x):
    xc = x*1.0
    xc.sort()
    sh = xc.shape[0]
    xmed = xc[sh//2]
    s0= int(sh*(1-.6827)/2)
    s1 = sh-s0
    x0 = xc[s0]
    x1 = xc[s1]
    s0b= int(sh*(1-.95)/2)
    s1b = sh-s0b
    try:
        x0b = xc[s0b]
        x1b = xc[s1b]
    except:
        x0b = xc[0]
        x1b = xc[-1]
    return((xmed,x0,x1,x0b,x1b))


# In[ ]:



# In[6]:


try:
    with open(folderstore+"gatinginfovslv31"+condition+shflag+".pickle","rb") as f:
        aps = pickle.load(f)
except:
    print('no final file!')
    try:
        with open(folderstore+"gatinginfovslv31_temp"+condition+shflag+".pickle","rb") as f:
            aps = pickle.load(f)
    except:
        print('no temporary file either!')
    
# with open("/mnt/data/gatinginfovslv2_temp.pickle","rb") as f:
#     aph = pickle.load(f)
    


# In[7]:


dsv = array([1,2,3,4,5,6,7,8,9])*7.5
dsvc = list(-1.0*dsv); dsvc.reverse()
dsv = concatenate((dsvc,[0],dsv))
print(dsv,dsv.shape)


# In[8]:


labels = ['Vspine','Vsoma','Vdendrite','Ca']
EL0s = [-70,-70,-70,0]
nr = 1
shle = {}
shle0 = {}
for lab in labels:
    shle[lab] = zeros((19,11))
    shle0[lab] = zeros((19,11))
    
#shle['Vspine'] = zeros((19,11))#shle['Vsoma'] = zeros((19,11))#shle['Vdendrite'] = zeros((19,11))#shle['Ca'] = zeros((19,11))

for il,lab in enumerate(labels):
    EL0s = list(aps[('EL0',0)])
    EL0 = EL0s[il]
    # columns: no inh, shaft inh, axo-spinal inh
    ap0 = aps[(lab,0.0,0)]-EL0
    for i in range(1,nr):
        ap0 = row_stack((ap0,aps[(lab,0.0,i)]-EL0))

    #Effect very close to inh
    r0shaft = ap0[:,1]/ap0[:,0]
    r0axspi = ap0[:,2]/ap0[:,0]
    for j,dss in enumerate(dsv):
        ap = aps[(lab,dss,0)]-EL0
        for i in range(1,nr):
            ap = row_stack((ap,aps[(lab,dss,i)]-EL0))
        rshaft = ap[:,1]/ap[:,0]
        raxspi = ap[:,2]/ap[:,0]
        difr = column_stack(((1-rshaft)/(1-r0shaft),(1-raxspi)/(1-r0axspi)))
        # difr = column_stack((rshaft,raxspi))
        shle0[lab][j,1:] = concatenate((getint(difr[:,0]),getint(difr[:,1])))
        shle0[lab][j,0] = dss
        difr[:,0] = difr[:,0]*(1-r0shaft)
        difr[:,1] = difr[:,1]*(1-r0axspi)
        shle[lab][j,1:] = concatenate((getint(difr[:,0]),getint(difr[:,1])))
        shle[lab][j,0] = dss


# In[11]:


shleh = shle
shleh0 = shle0


# In[18]:


with open(folderstore+"gatinginfo_lengthv31"+condition+shflag+".pickle","wb") as f:
    pickle.dump([shle,shleh,shle0,shleh0],f)
    


# In[ ]:




