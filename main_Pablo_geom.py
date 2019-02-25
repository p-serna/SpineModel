# 
#   Simulation Cable Eq. with compartments
#
#   Based on Smith et al. 2013 scripts.


import numpy as np
import pickle
import time
from numpy.random import exponential, randint
from numpy import *   #ones, cumsum, sum, isscalar
import sys

#from neuron import *


import PS_lib as lb
import PS_storage as st


if len(sys.argv)>1:
    il = int(sys.argv[1])
else:
    il = randint(0,100)
        
print("Estamos en dendrita",il,":")

#~ # Data is stored here      
data = st.dataStorage()


data.dt = 0.1
data.NMDA = False

lb.h.dt = data.dt
NMDA = data.NMDA

model = lb.loadNeuron("L23.hoc",axon=False)
for sec in lb.h.allsec():
    sec.L = sec.L*0.7
for sec in model.dend:
    sec.diam = sec.diam*0.7
    
#model.addDend(locus="dendA1",L=4.0,D=1.5)

model.addSpne(locus=il,ilocus=0.5,L=1.0,D=1.0,Lneck=1.4,Dneck=0.1)

#lb.h.topology()

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
vspneRec[1].record(sp(0.5)._ref_v)

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


nnmda = 8
nGABA = 35
gmaxN = 50*nnmda/1e6
gmaxG = nGABA*30e-6

lb.add_NMDAsyns(model, locs=[[0,0.5]], gmax=gmaxN,tau2=60.0)  
lb.add_GABAsyns(model, locs=[[0,1]], spne=True, gmax=gmaxG,tau1=0.8,tau2=20.0)  

NMDA = model.NMDAlist[0]
GABA = model.GABAlist[0]

NCN = model.ncNMDAlist[0]
NCG = model.ncGABAlist[0]

stimN=lb.h.NetStim();stimN.number = 1;
NCN = lb.h.NetCon(stimN,NMDA,0,0,gmaxN)
model.ncNMDAlist[0] = NCN

stimG=lb.h.NetStim();stimG.number = 1;
NCG = lb.h.NetCon(stimG,GABA,0,0,gmaxG)
model.ncGABAlist[0] = NCG

global model
lb.init_active(model,soma=True,dend=True,dendCa=True,spne=True,dendNa=True)




caDendRec = []
sp = model.spne[0]
caDendRec.append(lb.h.Vector())
caDendRec.append(lb.h.Vector())
caDendRec.append(lb.h.Vector())
caDendRec.append(lb.h.Vector())
caDendRec[0].record(sp(0.5)._ref_ica) 
#~ caDendRec[1].record(sp(0.5)._ref_it)
caDendRec[2].record(model.NMDAlist[0]._ref_i)
caDendRec[3].record(sp(1.0)._ref_cai) 

NC2 = model.ncNMDAlist[0]
NC2.delay = 0.
NCG = model.ncGABAlist[0]
NCG.delay = 3000


#~ lb.h.finitialize(model.E_PAS)
#~ lb.neuron.run(t_stop)
#~ plot(trec,vspneRec[0])
#~ plot(trec,vDendRec[0])
#~ figure(2)
#~ #plot(trec,caDendRec[0])
#~ plot(trec,caDendRec[3])

#%pylab

from matplotlib import *
from matplotlib.pyplot import *





tdels = 10**linspace(-0.5,log10(50),30)
tdels2 = hstack((-tdels[arange(len(tdels)-1,-1,-1)],0,tdels))

#-------------- VArying Neck resistance

NMDA.tau1 = 2.0
NMDA.tau2 = 40.0

spineArea =  sp.L*sp.diam+1.8*sp.diam**2/4 # um^2
CaTcond = 2000 # pS

NCG.weight[0] = 30*35/spineArea*1e-10*1e4
sp.gbar_itL = 5*CaTcond/spineArea*1e-10 # S/ cm^2   not pS/um^2 as in mod file
NCN.weight[0] = 50*0.0002/1e3


mes = []
mp = []
#Rnecks = 10**linspace(log10(1e6),log10(1e9),51)
Rnecks = 10**linspace(log10(1e6),log10(1e9),51)
Rnecks = 1e6*array([10,20,50,100,200,500,1000])
t_stop = 2000
ton = 1100
lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

for Rneck in Rnecks:
    neck = model.neck[0]
    #Rneck = neck.L*1e-6/(neck.diam*1e-6/2.0)**2/pi*neck.Ra/100
    neck.Ra = Rneck/neck.L*(neck.diam*1e-6/2.0)**2*pi*100/1e-6
    NC.delay = 3050
    NCN.delay = 3050
    NCG.delay = 3050.0
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    mCat =  caDendRec[3][-1]
    mV = vspneRec[0][-1]
    
    mp.append([Rneck,mV,mCat])
    for tdel in tdels2:
        NC.delay = ton-50
        NCN.delay = ton-50
        NCG.delay = ton+tdel-50.0

        lb.h.finitialize(model.E_PAS)
        lb.neuron.run(t_stop)
        mes.append([Rneck,tdel,max(vspneRec[0])-mV,max(caDendRec[3])-mCat])
    print(Rneck/1e6)


mes = array(mes)

savetxt("DwRneckGeom"+str(il).zfill(3)+".dat",mes)

cm1 = 'Blues'
cm2 = 'Oranges'


#--------------- Some plots: A mejorar
#figure(5)
xc = log10(array(Rnecks)); xc = xc-min(xc); xc = xc/max(xc)
cmapf = get_cmap(cm1)
cmapg = get_cmap(cm2)

figure(5)
mt = mes
for j,R2 in enumerate(Rnecks):
    sel = abs(mt[:,0]-R2)<1e-3
    plot(mt[sel,1],mt[sel,3]/max(mt[sel,3]),c=cmapf(xc[j]))
        
xlabel("Delta t")
ylabel("Calcium [ $\Delta$ F/F ]")
legend(title="R$_{neck}$ (M$\Omega$)")

    
twinx()
for j,R2 in enumerate(Rnecks):
    sel = abs(mt[:,0]-R2)<1e-3
    plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',c=cmapg(xc[j]))

ylabel("$\Delta$ V/V$_{max}$")

title("Dendrite :  "+str(il))

savefig("Geomdend_"+str(il).zfill(3)+".png")
