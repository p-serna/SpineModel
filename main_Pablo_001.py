# 
#   Simulation Cable Eq. with compartments
#
#   Based on Smith et al. 2013 scripts.
import numpy as np
import pickle
import time
from numpy.random import exponential, randint
from numpy import *   #ones, cumsum, sum, isscalar

#from neuron import *


import PS_lib as lb
import PS_storage as st



# Data is stored here      
data = st.dataStorage()


data.dt = 0.1
data.NMDA = False

lb.h.dt = data.dt
NMDA = data.NMDA

model = lb.loadNeuron(axon=False)
model.addDend(locus="dendA1",L=4.0,D=1.5)
model.addSpne(locus="Dend002",ilocus=0.5,L=1.0,D=1.0,Lneck=1.4,Dneck=0.1)

lb.h.topology()

trec, vrec = lb.h.Vector(), lb.h.Vector()

gRec, iRec, vDendRec, caDendRec, vSecRec, vspneRec = [], [], [], [], [], []
gNMDA_rec, iNMDA_rec = [], []
trec.record(lb.h._ref_t)
vrec.record(model.soma(0.5)._ref_v)

recDend =True
recSec = True
if recDend:
    n=0
    for dend in model.dend:
        vDendRec.append(lb.h.Vector())
        caDendRec.append(lb.h.Vector())
        vDendRec[n].record(dend(0.5)._ref_v)
        #caDendRec[n].record(dend(0.5)._ref_gna_na)
        n=1

if recSec:
    n = 0
    for sec in lb.h.allsec():
        for seg in sec.allseg():
            vSecRec.append(lb.h.Vector())
            vSecRec[n].record(seg._ref_v)
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

lb.h.celsius = model.temperature

vspneRec.append(lb.h.Vector())
vspneRec.append(lb.h.Vector())
sp = model.spne[0]
vspneRec[0].record(sp(0.5)._ref_v)
sp = model.neck[0]
vspneRec[0].record(sp(0.5)._ref_v)

t_stop = 200.0

#



model.AMPAlist = []
model.ncAMPAlist = []

AMPA = lb.h.Exp2Syn(1,sec = model.spne[0])

tau1  = 0.5
tau2 = 8.
AMPA.tau1 = tau1
AMPA.tau2 = tau2

nampa = 50
gmax = 15*nampa/1e6

stimE=lb.h.NetStim();stimE.number = 1; 

NC = lb.h.NetCon(stimE,AMPA,0,0,gmax)


    
model.AMPAlist.append(AMPA)
model.ncAMPAlist.append(NC)
NC.delay = 0

lb.h.finitialize(model.E_PAS)

lb.neuron.run(t_stop)

neck = model.neck[0]
Rneck = neck.L*1e-6/(neck.diam*1e-6/2.0)**2/pi*neck.Ra/100
dend = model.dend[0]
Rdend = dend.L*1e-6/2.0/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100
dend = model.dend[2]
Rdend +=dend.L*1e-6/2.0/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100

%pylab

plot(trec,vspneRec[0])
plot(trec,vDendRec[0])



nnmda = 8
gmax = 50*nnmda/1e6

global model
lb.init_active(model,soma=True,dend=True,dendCa=True,spne=True,dendNa=True)
lb.add_NMDAsyns(model, locs=[[0,1]], gmax=gmax,tau2=60.0)  




caDendRec = []
sp = model.spne[0]
caDendRec.append(lb.h.Vector())
caDendRec.append(lb.h.Vector())
caDendRec.append(lb.h.Vector())
caDendRec.append(lb.h.Vector())
caDendRec[0].record(sp(0.5)._ref_ica) 
caDendRec[1].record(sp(0.5)._ref_ik)
caDendRec[2].record(model.NMDAlist[0]._ref_i)
caDendRec[3].record(sp(1.0)._ref_cai) 

NC2 = model.ncNMDAlist[0]
NC2.delay = 0.


lb.h.finitialize(model.E_PAS)
lb.neuron.run(t_stop)
plot(trec,vspneRec[0])
plot(trec,vDendRec[0])
figure(2)
#plot(trec,caDendRec[0])
plot(trec,caDendRec[3])



mes = []

tdels = 10**linspace(-0.5,log10(50),30)
tdels2 = hstack((-tdels[arange(len(tdels)-1,-1,-1)],0,tdels))
t_stop = 2000
ton = 1100
for tdel in tdels2:
    lb.add_somaStim(model,p=0.5,onset=ton,dur=2, amp =1)
    NC.delay = ton+tdel-50
    NC2.delay = ton+tdel-50

    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    mes.append([tdel,max(vspneRec[0]),max(caDendRec[3])])
    print(tdel)

figure(1)
plot(trec,vspneRec[0],label="spine")
plot(trec,vDendRec[0],label="dend")
plot(trec,vrec,label="soma")
figure(2)
#plot(trec,caDendRec[0])
plot(trec,caDendRec[3])
