# 
#   Simulation Cable Eq. with compartments
#
#   Based on Smith et al. 2013 scripts.

import numpy as np
import pickle
import time
from numpy.random import exponential, randint
from numpy import *   #ones, cumsum, sum, isscalar
from matplotlib.pylab import * 
plotose = True

#from neuron import *


import PS_lib as lb
import PS_storage as st

#Logical Variables 
# Should we record in the dendrite?
recDend =True
# Should we record in ??
recSec = False
# Should we print stuff?
verbose = True

# Parameter definitions
# Data is stored here      
data = st.dataStorage() # some default parameters defined.
data.dt = 0.1
data.NMDA = False


# Definition of the model.
lb.h.dt = data.dt
NMDA = data.NMDA
model = lb.loadNeuron("Basic2.hoc",axon=False)
# Adding piece of dendritic branch with spine
model.addDend(name="DendE",locus="dendA1",L=4.0,D=1.5)
model.addSpne(locus="DendE",ilocus=0.5,L=1.0,D=1.0,Lneck=1.4,Dneck=0.1)
data.model = model.__dict__
# Temperature of the neuron
lb.h.celsius = model.temperature


# show topology
if verbose:
    print(lb.h.topology())


# Measurement Quantities
trec, vrec = lb.h.Vector(), lb.h.Vector()
gRec, iRec,  vspneRec = [], [], []
gNMDA_rec, iNMDA_rec = [], []
trec.record(lb.h._ref_t)
vrec.record(model.soma(0.5)._ref_v)


if recDend:
    #n=0
    vDendRec = []
    caDendRec = []
    #For all dendrites
    for dend in model.dend:
        #Adding vectors for Voltage, and Calcium
        vDendRec.append(lb.h.Vector())
        caDendRec.append(lb.h.Vector())
        # Placing recording at mid-point in the dendritic branch
        vDendRec[-1].record(dend(0.5)._ref_v)
        # NO CALCIUM!?!?!
#Probably better to organize them in a dictionary        

# Spine voltage recording stuff
vspneRec.append(lb.h.Vector())
vspneRec.append(lb.h.Vector())
sp = model.spne[0]
vspneRec[0].record(sp(0.5)._ref_v)
sp = model.neck[0]
vspneRec[1].record(sp(0.5)._ref_v)


if recSec:
    #n = 0
    vSecRec = []
    # This is for all segments?!
    for sec in lb.h.allsec():
        for seg in sec.allseg():
            vSecRec.append(lb.h.Vector())
            vSecRec[-1].record(seg._ref_v)
            n+=1


#~ if NMDA:        
    #~ for n in np.arange(0, len(model.NMDAlist)):
        #~ loc = model.NMDAlist[n].get_loc()
        #~ h.pop_section()                        
        #~ gNMDA_rec.append(h.Vector())
        #~ iNMDA_rec.append(h.Vector())
        #~ gNMDA_rec[n].record(model.NMDAlist[n]._ref_g)
        #~ iNMDA_rec[n].record(model.NMDAlist[n]._ref_i)
    #~ gRec.append(gNMDA_rec)
    #~ iRec.append(iNMDA_rec)



t_stop = 200.0

#

# Resistances calculation
neck = model.neck[0]
Rneck = neck.L*1e-6/(neck.diam*1e-6/2.0)**2/pi*neck.Ra/100
dend = model.dend[0]
Rdend = dend.L*1e-6/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100
dend = model.dend[1]
Rdend +=dend.L*1e-6/2.0/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100
dend = model.dend[-1]
Rdend +=dend.L*1e-6/2.0/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100


# Adding synapses:
#Excitatory synapse
model.AMPAlist = []
model.ncAMPAlist = []

AMPA = lb.h.Exp2Syn(1,sec = model.spne[0])
tau1  = 0.5
tau2 = 8.0
AMPA.tau1 = tau1
AMPA.tau2 = tau2

nampa = 50
gmax = 15*nampa/1e6
stimE=lb.h.NetStim();stimE.number = 1; 
NC = lb.h.NetCon(stimE,AMPA,0,0,gmax)

model.AMPAlist.append(AMPA)
model.ncAMPAlist.append(NC)
NC.delay = 10

# NMDA part
nnmda = 2
gmaxN = 50*nnmda/1e6
lb.add_NMDAsyns(model, locs=[[0,0.5]], gmax=gmaxN,tau2=20.0)  
NMDA = model.NMDAlist[0]
NCN = model.ncNMDAlist[0]
stimN=lb.h.NetStim();stimN.number = 1;
NCN = lb.h.NetCon(stimN,NMDA,0,0,gmaxN)
model.ncNMDAlist[0] = NCN


if plotose:
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)

    plot(trec,vspneRec[0])
    plot(trec,vDendRec[0])
    plot(trec,vrec)
    
ton = 50
nampas = arange(0,20)
mes = zeros((nampas.shape[0],5))
for i,nampa in enumerate(nampas):
    #lb.add_somaStim(model,p=0.5,onset=ton,dur=2, amp =0.0)
    NC.delay = ton-50
    NCN.delay = ton-50
    NC.weight[0] = nampa *15/1e6
    NCN.weight[0] = gmaxN

    # ~ print(NC.weight[0])
    lb.h.finitialize(model.E_PAS)
    # ~ print(NC.weight[0])

    lb.neuron.run(t_stop)
    
    #plot(trec,vspneRec[0])    
    
    current = abs((array(vDendRec[-2])-array(vDendRec[0]))/Rdend)

    mes[i,:] = [nampa,max(vspneRec[0]),max(vDendRec[1]),max(vrec),max(current)]

plot(mes[:,-1]/1e-12,mes[:,1])
plot(mes[:,-1]/1e-12,mes[:,2])
plot(mes[:,-1]/1e-12,mes[:,3])
xlabel("I(pA)")
