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

model = lb.loadNeuron("Basic2.hoc",axon=False)
model.addDend(name="DendE",locus="dendA1",L=4.0,D=1.5)
model.addSpne(locus="DendE",ilocus=0.5,L=1.0,D=1.0,Lneck=1.4,Dneck=0.1)

lb.h.topology()

trec, vrec = lb.h.Vector(), lb.h.Vector()

gRec, iRec, vDendRec, caDendRec, vSecRec, vspneRec = [], [], [], [], [], []
gNMDA_rec, iNMDA_rec = [], []
trec.record(lb.h._ref_t)
vrec.record(model.soma(0.5)._ref_v)

recDend =True
recSec = False
if recDend:
    n=0
    for dend in model.dend:
        vDendRec.append(lb.h.Vector())
        caDendRec.append(lb.h.Vector())
        vDendRec[n].record(dend(0.5)._ref_v)
        #caDendRec[n].record(dend(0.5)._ref_gna_na)
        n+=1

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

neck = model.neck[0]
Rneck = neck.L*1e-6/(neck.diam*1e-6/2.0)**2/pi*neck.Ra/100
dend = model.dend[0]
Rdend = dend.L*1e-6/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100
dend = model.dend[1]
Rdend +=dend.L*1e-6/2.0/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100
dend = model.dend[-1]
Rdend +=dend.L*1e-6/2.0/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100

%pylab


lb.h.finitialize(model.E_PAS)

lb.neuron.run(t_stop)

plot(trec,vspneRec[0])
plot(trec,vDendRec[0])


nnmda = 0.0003
nGABA = 35
gmaxN = 50*nnmda/1e3
gmaxG = nGABA*30e-3
lb.add_NMDAsyns(model, locs=[[0,0.5]], gmax=gmaxN,tau2=20.0)  
lb.add_GABAsyns(model, locs=[[0,1]], spne=True, gmax=gmaxG,tau1=1.5,tau2=20.0)  

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
caDendRec[1].record(model.NMDAlist[0]._ref_i)
caDendRec[2].record(model.NMDAlist[0]._ref_i)
caDendRec[3].record(sp(1.0)._ref_cai) 

NC2 = model.ncNMDAlist[0]
NC2.delay = 0.
NCG = model.ncGABAlist[0]
NCG.delay = 3000


lb.h.finitialize(model.E_PAS)
lb.neuron.run(t_stop)
plot(trec,vspneRec[0])
plot(trec,vDendRec[0])
figure(2)
#plot(trec,caDendRec[0])
plot(trec,caDendRec[3])


#sp.gbar_itL =  2000*1e-4
spineArea =  sp.L*sp.diam+1.8*sp.diam**2/4 # um^2
CaTcond = 2000 # pS
sp.gbar_itL = 20*CaTcond/spineArea*1e-10 # S/ cm^2   not pS/um^2 as in mod file

NMDA.tau1 = 2.0
NMDA.tau2 = 60.0
NCN.weight[0] = 50*0.0002/1e3
mes = []

tdels = 10**linspace(-0.5,log10(100),40)
tdels2 = hstack((-tdels[arange(len(tdels)-1,-1,-1)],0,tdels))
t_stop = 2000
ton = 1100
NCG.delay = 3000
for tdel in tdels2:
    lb.add_somaStim(model,p=0.5,onset=ton,dur=2, amp =1)
    NC.delay = ton-tdel-50
    NCN.delay = ton-tdel-50
    

    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    mes.append([tdel,max(vspneRec[0]),max(caDendRec[3])])
    print(tdel)
mes = array(mes)
messG = mes


figure(3)
plot(mes[:,0],mes[:,2])



ton = 3000
lb.add_somaStim(model,p=0.5,onset=ton,dur=2, amp =1)
NC.delay = 3050
NCN.delay = 3050
NCG.delay = 1000.0-50.
NCG.weight[0] = 30*35/spineArea*1e-10*1e4

lb.h.finitialize(model.E_PAS)
lb.neuron.run(t_stop)


figure(1)
plot(trec,vspneRec[0],label="spine")
plot(trec,vDendRec[0],label="dend")
plot(trec,vrec,label="soma")
figure(2)
#plot(trec,caDendRec[0])
plot(trec,caDendRec[3])


#--------------------------------------------------
#
#
#   Inh vs exc.

mes = []

NC.delay = 3050
NCN.delay = 3050
NCG.delay = 3050.0
lb.h.finitialize(model.E_PAS)
lb.neuron.run(t_stop)
mCat =  caDendRec[3][-1]
mV = vspneRec[0][-1]
    

tdels = 10**linspace(-0.5,log10(100),40)
tdels2 = hstack((-tdels[arange(len(tdels)-1,-1,-1)],0,tdels,200))
t_stop = 2000
ton = 1100

NCG.weight[0] = 30*35/spineArea*1e-10*1e4
sp.gbar_itL = 20*CaTcond/spineArea*1e-10 # S/ cm^2   not pS/um^2 as in mod file
NCN.weight[0] = 50*0.0002/1e3

lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)
for tdel in tdels2:
    NC.delay = ton-50
    NCN.delay = ton-50
    NCG.delay = ton+tdel-50.0

    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    mes.append([tdel,max(vspneRec[0])-mV,max(caDendRec[3])-mCat])
    print(tdel)
mes = array(mes)


figure(4)
plot(mes[:,0],(mes[:,2])/max(mes[:,2]),label="Calcium")
xlabel("Delta t")
ylabel("Ca$^{2+}_i$ [ ($\Delta$ F)/F ]")
#legend()
twinx()
plot(mes[:,0],(mes[:,1])/max(mes[:,1]),"C1--",label="Voltage")
ylabel("Voltage [$\Delta$ V/max($\Delta$ V)]")
sel = mes[:,1]/max(mes[:,1])==1
tp = mes[sel,0]
sel = tp>0 
vlines(tp[sel][0],0,1.2,'k',linestyles="dashed")
ylim(min(mes[:,1])/max(mes[:,1])*0.985,1.015)

#savefig("0T-Inhvsexc.png")

mes0 = mes


neck = model.neck[0]
Rneck = neck.L*1e-6/(neck.diam*1e-6/2.0)**2/pi*neck.Ra/100
dend = model.dend[0]
Rdend = dend.L*1e-6/2.0/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100
dend = model.dend[1]
Rdend +=dend.L*1e-6/2.0/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100
dend = model.dend[2]
Rdend +=dend.L*1e-6/2.0/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100

def cRdend():
    dend = model.dend[0]
    Rdend = dend.L*1e-6/2.0/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100
    dend = model.dend[1]
    Rdend +=dend.L*1e-6/2.0/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100
    dend = model.dend[2]
    Rdend +=dend.L*1e-6/2.0/(dend.diam*1e-6/2.0)**2/pi*dend.Ra/100
    return Rdend/1e6

tdels = 10**linspace(-0.5,log10(100),40)
tdels2 = hstack((-200,-tdels[arange(len(tdels)-1,-1,-1)],0,tdels,150,200))

#-------------- VArying Neck resistance

mes = []
mp = []
Rnecks = array([1,10,20,50,100,200,500,1000,10000])*1e6

t_stop = 2000
ton = 1100

# There is AP!
lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

NCG.weight[0] = 30*35/spineArea*1e-10*1e4
sp.gbar_itL = 20*CaTcond/spineArea*1e-10 # S/ cm^2   not pS/um^2 as in mod file
NCN.weight[0] = 50*0.0002/1e3

for Rneck in Rnecks:
    neck = model.neck[0]
    #Rneck = neck.L*1e-6/(neck.diam*1e-6/2.0)**2/pi*neck.Ra/100
    neck.Ra = Rneck/neck.L*(neck.diam*1e-6/2.0)**2*pi*100/1e-6
    lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

    NC.delay = 3050
    NCN.delay = 3050
    NCG.delay = 3050.0
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    mCat =  caDendRec[3][-1]
    mV = vspneRec[0][-1]
    
    lb.add_somaStim(model,p=0.5,onset=ton-500,dur=2, amp =2)
    
    NC.delay = 1050
    NCN.delay = 1050
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    maCat = max(caDendRec[3])-mCat
    
    mp.append([Rneck,mV,mCat])
    for tdel in tdels2:
        lb.add_somaStim(model,p=0.5,onset=ton+tdel,dur=2, amp =2)
        NC.delay = ton-50
        NCN.delay = ton-50
        NCG.delay = ton+1e5
        vca = []
        ca = caDendRec[3]
        for v in ca: vca.append((v-mCat)/maCat)
        vca = array(vca)
        t1 = sum(vca>1.01)*0.1
        t2 = sum(vca>1.2)*0.1
        t3 = sum(vca>1.4)*0.1
        
        lb.h.finitialize(model.E_PAS)
        lb.neuron.run(t_stop)
        if t1>0:
            mes.append([Rneck,tdel,t1,max(caDendRec[3])-mCat,t2/t1,t3/t1])
        else:
            mes.append([Rneck,tdel,t1,max(caDendRec[3])-mCat,0.,0.])
    print(Rneck/1e6)


mes = array(mes)
savetxt("DbAPRneckRdend36MOhm.dat",mes)

#--------- Paleta de colores
cm = 'Blues'
cm2 = 'Oranges'
xc = log10(array(Rnecks)); xc = xc-min(xc); xc = xc/max(xc)
cmapf = get_cmap(cm)
cmapg = get_cmap(cm2)

figure(5)

mt = mes
for j,R2 in enumerate(Rnecks[1:-1]):
    sel = abs(mt[:,0]-R2)<1e-3
    plot(mt[sel,1],mt[sel,3]/mt[sel,3][-1],label=R2/1e6,c=cmapf(xc[j]))
        
xlabel("$\Delta$ t")
ylabel("Calcium [ $\Delta$ F/F ]")
    #legend()    
#legend(title="R$_{neck}$ (M$\Omega$)")
ylim(0,1.7)
twinx()
for j,R2 in enumerate(Rnecks[1:-1]):
    sel = abs(mt[:,0]-R2)<1e-3
    plot(mt[sel,1],-mt[sel,2]*(1-(mt[sel,5]>0))/max(mt[sel,2])*0.3+mt[sel,5]*0.6,'--',label=R2/1e6,c=cmapg(xc[j]))
    #plot(mes[:,1],-mes[:,2]*(1-(mes[:,4]>0))/max(mes[:,2])+mes[:,4]*0.5)
    #plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',label=R2/1e6,c=cmapg(xc[j]))

ylabel("$\Delta$ w")

title("Z$_{dend}$ = "+str(int(Rdend/1e6))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
ylim(-0.3,0.4)
savefig("STDPzdend"+str(int(Rdend/1e6)).zfill(4)+"mOhm.png")

close(5)
figure(5)
for j,R2 in enumerate(Rnecks[1:-1]):
    sel = abs(mt[:,0]-R2)<1e-3
    plot(mt[sel,1],-mt[sel,2]*(1-(mt[sel,5]>0))/max(mt[sel,2])*0.3+mt[sel,5]*0.6,'--',label=R2/1e6,c=cmapg(xc[j]))
    #plot(mes[:,1],-mes[:,2]*(1-(mes[:,4]>0))/max(mes[:,2])+mes[:,4]*0.5)
    #plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',label=R2/1e6,c=cmapg(xc[j]))

ylabel("$\Delta$ w$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(int(Rdend/1e6))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
ylim(-0.3,0.4)
savefig("STDPzdend"+str(int(Rdend/1e6)).zfill(4)+"mOhmA.png")



# Varying both neck and dendritic resistance.

mes = []
mp = []
#Rnecks = 10**linspace(log10(1e6),log10(1e9),21)

Rnecks = array([1,10,20,50,100,200,500,1000,2000,10000])*1e6
Rdends0 = array([70,100,200,300,400,500,750,1000,1500,2000])*1e6


Rnecks = array([50,500])*1e6
Rdends0 = array([70,200,500])*1e6
Rs = [[i,j] for i in Rnecks for j in Rdends0]

t_stop = 2000
ton = 1100
lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

dendt = model.dend[1]
dendt.L = 50.0
Rmin = 70.84852670568594e6 #177.12856043961426

Rdends=[]
for Rdend in Rdends0:
    dend = model.dend[0]
    dend.Ra = max(Rdend-Rmin,1e-1)/(dend.L/2.0)*(dend.diam*1e-6/2.0)**2*pi*100/1e-6
    #print(Rdend-Rmin,dend.Ra)
    Rdend2 = cRdend()
    Rdends.append(Rdend2)

dendt = model.dend[1]
for Rneck,Rdend in Rs:
    neck = model.neck[0]
    dend = model.dend[0]
    #Rneck = neck.L*1e-6/(neck.diam*1e-6/2.0)**2/pi*neck.Ra/100
    neck.Ra = Rneck/neck.L*(neck.diam*1e-6/2.0)**2*pi*100/1e-6
    dend.Ra = max(Rdend-Rmin,1e-1)/(dend.L/2.0)*(dend.diam*1e-6/2.0)**2*pi*100/1e-6
    
    Rdend2 = cRdend()
    #print(Rdend2)
    lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

    NC.delay = 3050
    NCN.delay = 3050
    NCG.delay = 3050.0
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    mCat =  caDendRec[3][-1]
    mV = vspneRec[0][-1]
    
    lb.add_somaStim(model,p=0.5,onset=ton-500,dur=2, amp =2)
    
    NC.delay = 1050
    NCN.delay = 1050
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    maCat = max(caDendRec[3])-mCat
    
    mp.append([Rdend,Rneck,mV,mCat])
    for tdel in tdels2:
        lb.add_somaStim(model,p=0.5,onset=ton+tdel,dur=2, amp =2)
        NC.delay = ton-50
        NCN.delay = ton-50
        NCG.delay = ton+1e5
        
        lb.h.finitialize(model.E_PAS)
        lb.neuron.run(t_stop)
        ca = caDendRec[3]
        vca = []
        for v in ca: vca.append((v-mCat)/maCat)
        vca = array(vca)
        t1 = sum(vca>1.01)*0.1
        t2 = sum(vca>1.2)*0.1
        t3 = sum(vca>1.4)*0.1
        if t1>0:
            mes.append([Rdend2,Rneck,tdel,t1,max(caDendRec[3])-mCat,t2/t1,t3/t1])
        else:
            mes.append([Rdend2,Rneck,tdel,t1,max(caDendRec[3])-mCat,0.,0.0])
    print(Rdend/1e6,Rneck/1e6)


mes = array(mes)

mesR2s = mes
savetxt("DbAPwRR.dat",mes)
#
#savetxt("DbAPwR2s.dat",mesR2s)

cm = 'Blues'
cm2 = 'Oranges'


#--------------- Some plots: A mejorar

figure(5)
xc = log10(array(Rnecks)); xc = xc-min(xc); xc = xc/max(xc)
cmapf = get_cmap(cm)
cmapg = get_cmap(cm2)

#ies = [0,5,10,15,17,20]
ies = arange(len(Rdends))
#Rdends = [mes[1,0],mes[81,0],mes[162,0],mes[81*3,0],mes[81*4,0],mes[81*5,0],mes[81*6,0],mes[81*7,0],mes[81*8,0]]

for it in ies:
    R = Rdends[it]
    figure(it)
    sel1 = abs(mes[:,0]-R)<1e-3
    mt = mes[sel1,:]
    for j,R2 in enumerate(Rnecks[1:-1]):
        sel = abs(mt[:,1]-R2)<1e-3
        plot(mt[sel,2],mt[sel,4]/max(mt[sel,4]),label=R2/1e6,c=cmapf(xc[j]))
            
    xlabel("Delta t")
    ylabel("Calcium [ $\Delta$ F/F ]")
    legend(title="R$_{Neck}$ (M$\Omega$)")    
        
    twinx()
    for j,R2 in enumerate(Rnecks[1:-1]):
        sel = abs(mt[:,1]-R2)<1e-3
        plot(mt[sel,2],mt[sel,3]/max(mt[sel,3]),'--',label=R2/1e6,c=cmapg(xc[j]))

    ylabel("$\Delta$ V/V$_{max}$")

    title("Z$_{dend}$ = "+str(int(R))+" M $\Omega$")

    savefig("zdend"+str(int(R)).zfill(4)+"mOhm.png")

figure(5)
#~ xc = log10(array(Rnecks)); xc = xc-min(xc); xc = xc/max(xc)
#~ cmapf = get_cmap(cm)
#~ cmapg = get_cmap(cm2)

#ies = [0,5,10,15,17,20]
ies = arange(len(Rnecks))
for it in ies:
    R = Rnecks[it]
    #close(5)
    figure(it)
    sel1 = abs(mes[:,1]-R)<1e-3
    mt = mes[sel1,:]
    for j,R2 in enumerate(Rdends):
        sel = abs(mt[:,0]-R2)<1e-3
        plot(mt[sel,2],mt[sel,4]/max(mt[sel,4]),label=round(R2),c=cmapf(xc[j]))
            
    xlabel("Delta t")
    ylabel("Calcium [ $\Delta$ F/F ]")
        #legend()    
    legend(title="R$_{Dend}$ (M$\Omega$)")    
     
    twinx()
    for j,R2 in enumerate(Rdends):
        sel = abs(mt[:,0]-R2)<1e-3
        plot(mt[sel,2],mt[sel,3]/max(mt[sel,3]),'--',label=round(R2),c=cmapg(xc[j]))

    ylabel("$\Delta$ V/V$_{max}$")

    title("R$_{neck}$ = "+str(int(R/1e6))+" M $\Omega$")

    savefig("Rneck"+str(int(R/1e6)).zfill(4)+"mOhm.png")


#~ #ies = [0,5,10,15,17,20]
#~ for it in ies:
    #~ R = Rnecks[it]
    #~ close(5)
    #~ figure(5)
    #~ sel1 = abs(mes[:,1]-R)<1e-3
    #~ mt = mes[sel1,:]
    #~ for j,R2 in enumerate(Rnecks):
        #~ sel = abs(mt[:,0]-R2)<1e-3
        #~ plot(mt[sel,2],mt[sel,4]/max(mt[sel,4]),label=R/1e6,c=cmapf(xc[j]))
            
    #~ xlabel("Delta t")
    #~ ylabel("Calcium [ $\Delta$ F/F ]")
        #~ #legend()    
        
    #~ twinx()
    #~ for j,R2 in enumerate(Rnecks):
        #~ sel = abs(mt[:,0]-R2)<1e-3
        #~ plot(mt[sel,2],mt[sel,3]/max(mt[sel,3]),'--',label=R/1e6,c=cmapg(xc[j]))

    #~ ylabel("$\Delta$ V/V$_{max}$")

#legend()
#~ figure(1)
#~ plot(trec,vspneRec[0],label="spine")
#~ plot(trec,vDendRec[0],label="dend")
#~ plot(trec,vrec,label="soma")
#~ figure(2)
#~ #plot(trec,caDendRec[0])
#~ plot(trec,caDendRec[3])



#------------- with GABA





# Varying both neck and dendritic resistance.

mes = []
mp = []
#Rnecks = 10**linspace(log10(1e6),log10(1e9),21)

Rnecks = array([1,10,20,50,100,200,500,1000,2000,10000])*1e6
Rdends0 = array([70,100,200,300,400,500,750,1000,1500,2000])*1e6

Rnecks = array([50,500])*1e6
Rdends0 = array([70,200,500])
Rs = [[i,j] for i in Rnecks for j in Rdends0]

t_stop = 2000
ton = 1100
lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

dendt = model.dend[1]
dendt.L = 50.0
Rmin = 70.84852670568594e6 #177.12856043961426

Rdends=[]
for Rdend in Rdends0:
    dend = model.dend[0]
    dend.Ra = max(Rdend-Rmin,1e-1)/(dend.L/2.0)*(dend.diam*1e-6/2.0)**2*pi*100/1e-6
    #print(Rdend-Rmin,dend.Ra)
    Rdend2 = cRdend()
    Rdends.append(Rdend2)

dendt = model.dend[1]
for Rneck,Rdend in Rs:
    neck = model.neck[0]
    dend = model.dend[0]
    #Rneck = neck.L*1e-6/(neck.diam*1e-6/2.0)**2/pi*neck.Ra/100
    neck.Ra = Rneck/neck.L*(neck.diam*1e-6/2.0)**2*pi*100/1e-6
    dend.Ra = max(Rdend-Rmin,1e-1)/(dend.L/2.0)*(dend.diam*1e-6/2.0)**2*pi*100/1e-6
    
    Rdend2 = cRdend()
    #print(Rdend2)
    lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

    NC.delay = 3050
    NCN.delay = 3050
    NCG.delay = 3050.0
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    mCat =  caDendRec[3][-1]
    mV = vspneRec[0][-1]
    
    lb.add_somaStim(model,p=0.5,onset=ton-500,dur=2, amp =2)
    
    NC.delay = 1050
    NCN.delay = 1050
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    maCat = max(caDendRec[3])-mCat
    
    mp.append([Rdend,Rneck,mV,mCat])
    for tdel in tdels2:
        lb.add_somaStim(model,p=0.5,onset=ton+tdel,dur=2, amp =2)
        NC.delay = ton-50
        NCN.delay = ton-50
        NCG.delay = ton+tdel-50.0-5.0
        
        
        lb.h.finitialize(model.E_PAS)
        lb.neuron.run(t_stop)
        
        vca = []
        ca = caDendRec[3]
        for v in ca: vca.append((v-mCat)/maCat)
        vca = array(vca)
        t1 = sum(vca>1.01)*0.1
        t2 = sum(vca>1.2)*0.1
        t3 = sum(vca>1.4)*0.1
        if t1>0:
            mes.append([Rdend2,Rneck,tdel,t1,max(caDendRec[3])-mCat,t2/t1,t3/t1])
        else:
            mes.append([Rdend2,Rneck,tdel,t1,max(caDendRec[3])-mCat,0.,0.0])
    print(Rdend/1e6,Rneck/1e6)


mes = array(mes)

savetxt("DbAPwR2sGm5.dat",mes)

mesR2sGm5 = mes




# Varying both neck and dendritic resistance.

mes = []
mp = []
#Rnecks = 10**linspace(log10(1e6),log10(1e9),21)

Rnecks = array([1,10,20,50,100,200,500,1000,2000,10000])*1e6
Rdends0 = array([70,100,200,300,400,500,750,1000,1500,2000])*1e6
Rs = [[i,j] for i in Rnecks for j in Rdends0]

t_stop = 2000
ton = 1100
lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

dendt = model.dend[1]
dendt.L = 50.0
Rmin = 70.84852670568594e6 #177.12856043961426

Rdends=[]
for Rdend in Rdends0:
    dend = model.dend[0]
    dend.Ra = max(Rdend-Rmin,1e-1)/(dend.L/2.0)*(dend.diam*1e-6/2.0)**2*pi*100/1e-6
    #print(Rdend-Rmin,dend.Ra)
    Rdend2 = cRdend()
    Rdends.append(Rdend2)

dendt = model.dend[1]
for Rneck,Rdend in Rs:
    neck = model.neck[0]
    dend = model.dend[0]
    #Rneck = neck.L*1e-6/(neck.diam*1e-6/2.0)**2/pi*neck.Ra/100
    neck.Ra = Rneck/neck.L*(neck.diam*1e-6/2.0)**2*pi*100/1e-6
    dend.Ra = max(Rdend-Rmin,1e-1)/(dend.L/2.0)*(dend.diam*1e-6/2.0)**2*pi*100/1e-6
    
    Rdend2 = cRdend()
    #print(Rdend2)
    lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

    NC.delay = 3050
    NCN.delay = 3050
    NCG.delay = 3050.0
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    mCat =  caDendRec[3][-1]
    mV = vspneRec[0][-1]
    
    lb.add_somaStim(model,p=0.5,onset=ton-500,dur=2, amp =2)
    
    NC.delay = 1050
    NCN.delay = 1050
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    maCat = max(caDendRec[3])-mCat
    
    mp.append([Rdend,Rneck,mV,mCat])
    for tdel in tdels2:
        lb.add_somaStim(model,p=0.5,onset=ton+tdel,dur=2, amp =2)
        NC.delay = ton-50
        NCN.delay = ton-50
        NCG.delay = ton+tdel-50.0+5.0
        vca = []
        ca = caDendRec[3]
        for v in ca: vca.append((v-mCat)/maCat)
        vca = array(vca)
        t1 = sum(vca>1.01)*0.1
        t2 = sum(vca>1.2)*0.1
        t3 = sum(vca>1.4)*0.1
        
        lb.h.finitialize(model.E_PAS)
        lb.neuron.run(t_stop)
        if t1>0:
            mes.append([Rdend2,Rneck,tdel,t1,max(caDendRec[3])-mCat,t2/t1,t3/t1])
        else:
            mes.append([Rdend2,Rneck,tdel,t1,max(caDendRec[3])-mCat,0.,0.0])
    print(Rdend/1e6,Rneck/1e6)


mes = array(mes)

savetxt("DbAPwR2sGp5.dat",mes)

mesR2sGp5 = mes




# Varying both neck and dendritic resistance.

mes = []
mp = []
#Rnecks = 10**linspace(log10(1e6),log10(1e9),21)

Rnecks = array([1,10,20,50,100,200,500,1000,2000,10000])*1e6
Rdends0 = array([70,100,200,300,400,500,750,1000,1500,2000])*1e6
Rs = [[i,j] for i in Rnecks for j in Rdends0]

t_stop = 2000
ton = 1100
lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

dendt = model.dend[1]
dendt.L = 50.0
Rmin = 70.84852670568594e6 #177.12856043961426

Rdends=[]
for Rdend in Rdends0:
    dend = model.dend[0]
    dend.Ra = max(Rdend-Rmin,1e-1)/(dend.L/2.0)*(dend.diam*1e-6/2.0)**2*pi*100/1e-6
    #print(Rdend-Rmin,dend.Ra)
    Rdend2 = cRdend()
    Rdends.append(Rdend2)

dendt = model.dend[1]
for Rneck,Rdend in Rs:
    neck = model.neck[0]
    dend = model.dend[0]
    #Rneck = neck.L*1e-6/(neck.diam*1e-6/2.0)**2/pi*neck.Ra/100
    neck.Ra = Rneck/neck.L*(neck.diam*1e-6/2.0)**2*pi*100/1e-6
    dend.Ra = max(Rdend-Rmin,1e-1)/(dend.L/2.0)*(dend.diam*1e-6/2.0)**2*pi*100/1e-6
    
    Rdend2 = cRdend()
    #print(Rdend2)
    lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

    NC.delay = 3050
    NCN.delay = 3050
    NCG.delay = 3050.0
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    mCat =  caDendRec[3][-1]
    mV = vspneRec[0][-1]
    
    lb.add_somaStim(model,p=0.5,onset=ton-500,dur=2, amp =2)
    
    NC.delay = 1050
    NCN.delay = 1050
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    maCat = max(caDendRec[3])-mCat
    
    mp.append([Rdend,Rneck,mV,mCat])
    for tdel in tdels2:
        lb.add_somaStim(model,p=0.5,onset=ton+tdel,dur=2, amp =2)
        NC.delay = ton-50
        NCN.delay = ton-50
        NCG.delay = ton-50.0-5.0
        vca = []
        ca = caDendRec[3]
        for v in ca: vca.append((v-mCat)/maCat)
        vca = array(vca)
        t1 = sum(vca>1.01)*0.1
        t2 = sum(vca>1.2)*0.1
        t3 = sum(vca>1.4)*0.1
        
        lb.h.finitialize(model.E_PAS)
        lb.neuron.run(t_stop)
        if t1>0:
            mes.append([Rdend2,Rneck,tdel,t1,max(caDendRec[3])-mCat,t2/t1,t3/t1])
        else:
            mes.append([Rdend2,Rneck,tdel,t1,max(caDendRec[3])-mCat,0.,0.0])
    print(Rdend/1e6,Rneck/1e6)


mes = array(mes)

savetxt("DbAPwR2sGprem5.dat",mes)

mesR2sGprem5 = mes


# Varying both neck and dendritic resistance.

mes = []
mp = []
#Rnecks = 10**linspace(log10(1e6),log10(1e9),21)

Rnecks = array([1,10,20,50,100,200,500,1000,2000,10000])*1e6
Rdends0 = array([70,100,200,300,400,500,750,1000,1500,2000])*1e6
Rs = [[i,j] for i in Rnecks for j in Rdends0]

t_stop = 2000
ton = 1100
lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

dendt = model.dend[1]
dendt.L = 50.0
Rmin = 70.84852670568594e6 #177.12856043961426

Rdends=[]
for Rdend in Rdends0:
    dend = model.dend[0]
    dend.Ra = max(Rdend-Rmin,1e-1)/(dend.L/2.0)*(dend.diam*1e-6/2.0)**2*pi*100/1e-6
    #print(Rdend-Rmin,dend.Ra)
    Rdend2 = cRdend()
    Rdends.append(Rdend2)

dendt = model.dend[1]
for Rneck,Rdend in Rs:
    neck = model.neck[0]
    dend = model.dend[0]
    #Rneck = neck.L*1e-6/(neck.diam*1e-6/2.0)**2/pi*neck.Ra/100
    neck.Ra = Rneck/neck.L*(neck.diam*1e-6/2.0)**2*pi*100/1e-6
    dend.Ra = max(Rdend-Rmin,1e-1)/(dend.L/2.0)*(dend.diam*1e-6/2.0)**2*pi*100/1e-6
    
    Rdend2 = cRdend()
    #print(Rdend2)
    lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

    NC.delay = 3050
    NCN.delay = 3050
    NCG.delay = 3050.0
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    mCat =  caDendRec[3][-1]
    mV = vspneRec[0][-1]
    
    lb.add_somaStim(model,p=0.5,onset=ton-500,dur=2, amp =2)
    
    NC.delay = 1050
    NCN.delay = 1050
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    maCat = max(caDendRec[3])-mCat
    
    mp.append([Rdend,Rneck,mV,mCat])
    for tdel in tdels2:
        lb.add_somaStim(model,p=0.5,onset=ton+tdel,dur=2, amp =2)
        NC.delay = ton-50
        NCN.delay = ton-50
        NCG.delay = ton-50.0+5.0
        vca = []
        ca = caDendRec[3]
        for v in ca: vca.append((v-mCat)/maCat)
        vca = array(vca)
        t1 = sum(vca>1.01)*0.1
        t2 = sum(vca>1.2)*0.1
        t3 = sum(vca>1.4)*0.1
        
        lb.h.finitialize(model.E_PAS)
        lb.neuron.run(t_stop)
        if t1>0:
            mes.append([Rdend2,Rneck,tdel,t1,max(caDendRec[3])-mCat,t2/t1,t3/t1])
        else:
            mes.append([Rdend2,Rneck,tdel,t1,max(caDendRec[3])-mCat,0.,0.0])
    print(Rdend/1e6,Rneck/1e6)


mes = array(mes)

savetxt("DbAPwR2sGprep5.dat",mes)

mesR2sGprep5 = mes


#---------------------------------------------------------------
#
#
#


for i in range(10): close()


it = 2
xc = log10(array(Rnecks[1:-2])); xc = xc-min(xc); xc = xc/max(xc)

mes = mesR2sGm5
R = Rdends[it]
figure(it)
sel1 = abs(mes[:,0]-R)<1e-3
mt = mes[sel1,:]

figure(5)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    plot(mt[sel,2],mt[sel,4]/mt[sel,4][0],label=R2/1e6,c=cmapf(xc[j]))
    #plot(mt[sel,2],-mt[sel,3]*(1-(mt[sel,6]>0))/max(mt[sel,3])*0.3+mt[sel,6]*0.6,'--',label=R2/1e6,c=cmapg(xc[j]))
    #plot(mes[:,1],-mes[:,2]*(1-(mes[:,4]>0))/max(mes[:,2])+mes[:,4]*0.5)
    #plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',label=R2/1e6,c=cmapg(xc[j]))

ylabel("$\Delta$ w")
xlabel("$\Delta t$")
title("Z$_{dend}$ = "+str(int(Rdend/1e6))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
#ylim(-0.3,0.4)
savefig("0T-STDP_CA_GABAposm5.png")

#mes = mesR2sGprep5

figure(10)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    mtp = mt[sel,4]/mt[sel,4][0]
    plot(mt[sel,2],(mtp-1.25)*(mtp>1.25)*2-(mtp-1.1)*(mtp>1.1),label=R2/1e6,c=cmapf(xc[j]))

ylabel("$\Delta w$")
xlabel("$\Delta$ t")
title("Z$_{dend}$ = "+str(int(Rdend/1e6))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
ylim(-0.2,0.4)
xlim(-100,100)
savefig("0T-STDP_GABAposm5.png")


for i in range(10): close()



it = 2
xc = log10(array(Rnecks[1:-2])); xc = xc-min(xc); xc = xc/max(xc)

mes = mesR2sGp5
R = Rdends[it]
#figure(it)
sel1 = abs(mes[:,0]-R)<1e-3
mt = mes[sel1,:]

figure(5)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    plot(mt[sel,2],mt[sel,4]/mt[sel,4][0],label=R2/1e6,c=cmapf(xc[j]))
    #plot(mt[sel,2],-mt[sel,3]*(1-(mt[sel,6]>0))/max(mt[sel,3])*0.3+mt[sel,6]*0.6,'--',label=R2/1e6,c=cmapg(xc[j]))
    #plot(mes[:,1],-mes[:,2]*(1-(mes[:,4]>0))/max(mes[:,2])+mes[:,4]*0.5)
    #plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',label=R2/1e6,c=cmapg(xc[j]))

ylabel("$\Delta w$")
xlabel("$\Delta$ t")
title("Z$_{dend}$ = "+str(int(Rdend/1e6))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
#ylim(-0.3,0.4)
savefig("0T-STDP_CA_GABAposp5.png")

#mes = mesR2sGprep5

figure(10)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    mtp = mt[sel,4]/mt[sel,4][0]
    plot(mt[sel,2],(mtp-1.25)*(mtp>1.25)*2-(mtp-1.1)*(mtp>1.1),label=R2/1e6,c=cmapf(xc[j]))

ylabel("$\Delta$ w$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(int(Rdend/1e6))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
ylim(-0.2,0.4)
xlim(-100,100)
savefig("0T-STDP_GABAposp5.png")



for i in range(10): close()




it = 2
xc = log10(array(Rnecks[1:-2])); xc = xc-min(xc); xc = xc/max(xc)

mes = mesR2sGprep5
R = Rdends[it]
#figure(it)
sel1 = abs(mes[:,0]-R)<1e-3
mt = mes[sel1,:]

figure(5)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    plot(mt[sel,2],mt[sel,4]/mt[sel,4][0],label=R2/1e6,c=cmapf(xc[j]))
    #plot(mt[sel,2],-mt[sel,3]*(1-(mt[sel,6]>0))/max(mt[sel,3])*0.3+mt[sel,6]*0.6,'--',label=R2/1e6,c=cmapg(xc[j]))
    #plot(mes[:,1],-mes[:,2]*(1-(mes[:,4]>0))/max(mes[:,2])+mes[:,4]*0.5)
    #plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',label=R2/1e6,c=cmapg(xc[j]))

ylabel("$\Delta$ w$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(int(Rdend/1e6))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
#ylim(-0.3,0.4)
savefig("0T-STDP_CA_GABAprep5.png")

#mes = mesR2sGprep5

figure(10)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    mtp = mt[sel,4]/mt[sel,4][0]
    plot(mt[sel,2],(mtp-1.25)*(mtp>1.25)*2-(mtp-1.1)*(mtp>1.1),label=R2/1e6,c=cmapf(xc[j]))

ylabel("$\Delta$ w$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(int(Rdend/1e6))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
ylim(-0.2,0.4)
xlim(-100,100)
savefig("0T-STDP_GABAprep5.png")




for i in range(10): close()




it = 2
xc = log10(array(Rnecks[1:-2])); xc = xc-min(xc); xc = xc/max(xc)

mes = mesR2sGprem5
R = Rdends[it]
#figure(it)
sel1 = abs(mes[:,0]-R)<1e-3
mt = mes[sel1,:]

figure(5)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    plot(mt[sel,2],mt[sel,4]/mt[sel,4][0],label=R2/1e6,c=cmapf(xc[j]))
    #plot(mt[sel,2],-mt[sel,3]*(1-(mt[sel,6]>0))/max(mt[sel,3])*0.3+mt[sel,6]*0.6,'--',label=R2/1e6,c=cmapg(xc[j]))
    #plot(mes[:,1],-mes[:,2]*(1-(mes[:,4]>0))/max(mes[:,2])+mes[:,4]*0.5)
    #plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',label=R2/1e6,c=cmapg(xc[j]))

ylabel("$\Delta$ w$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(int(Rdend/1e6))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
#ylim(-0.3,0.4)
savefig("0T-STDP_CA_GABAprem5.png")

#mes = mesR2sGprep5

figure(10)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    mtp = mt[sel,4]/mt[sel,4][0]
    plot(mt[sel,2],(mtp-1.25)*(mtp>1.25)*2-(mtp-1.1)*(mtp>1.1),label=R2/1e6,c=cmapf(xc[j]))

ylabel("$\Delta$ w$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(int(Rdend/1e6))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
ylim(-0.2,0.4)
xlim(-100,100)
savefig("0T-STDP_GABAprem5.png")




RnecksA = Rnecks

Rnecks = Rnecks[[2,3,6,8,9]]

it = 2
xc = log10(array(Rnecks[1:-2])); xc = xc-min(xc); xc = xc/max(xc)

mes = mesR2sGm5
R = Rdends[it]
figure(it)
sel1 = abs(mes[:,0]-R)<1e-3
mt = mes[sel1,:]

figure(5)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    plot(mt[sel,2],mt[sel,4]/mt[sel,4][0],label=R2/1e6)
    #plot(mt[sel,2],-mt[sel,3]*(1-(mt[sel,6]>0))/max(mt[sel,3])*0.3+mt[sel,6]*0.6,'--',label=R2/1e6,c=cmapg(xc[j]))
    #plot(mes[:,1],-mes[:,2]*(1-(mes[:,4]>0))/max(mes[:,2])+mes[:,4]*0.5)
    #plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',label=R2/1e6,c=cmapg(xc[j]))

ylabel("$\Delta$ F/F$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(round(R))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
#ylim(-0.3,0.4)
savefig("1T-STDP_CA_GABAposm5.png")

#mes = mesR2sGprep5

figure(10)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    mtp = mt[sel,4]/mt[sel,4][0]
    plot(mt[sel,2],(mtp-1.25)*(mtp>1.25)*2-(mtp-1.1)*(mtp>1.1),label=R2/1e6)
mtt = mt[sel,:]
sel = array([10,63])
mtp = mtp[sel]
mtt = mtt[sel,2]
plot(mtt,(mtp-1.25)*(mtp>1.25)*2-(mtp-1.1)*(mtp>1.1),'C1:')
ylabel("$\Delta$ w$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(round(R))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
ylim(-0.2,0.45)
xlim(-100,100)
savefig("1T-STDP_GABAposm5.png")


for i in range(10): close()



it = 2
xc = log10(array(Rnecks[1:-2])); xc = xc-min(xc); xc = xc/max(xc)

mes = mesR2sGp5
R = Rdends[it]
#figure(it)
sel1 = abs(mes[:,0]-R)<1e-3
mt = mes[sel1,:]

figure(5)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    plot(mt[sel,2],mt[sel,4]/mt[sel,4][0],label=R2/1e6)
    #plot(mt[sel,2],-mt[sel,3]*(1-(mt[sel,6]>0))/max(mt[sel,3])*0.3+mt[sel,6]*0.6,'--',label=R2/1e6,c=cmapg(xc[j]))
    #plot(mes[:,1],-mes[:,2]*(1-(mes[:,4]>0))/max(mes[:,2])+mes[:,4]*0.5)
    #plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',label=R2/1e6,c=cmapg(xc[j]))

ylabel("$\Delta$ F/F$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(round(R))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
#ylim(-0.3,0.4)
savefig("1T-STDP_CA_GABAposp5.png")

#mes = mesR2sGprep5

figure(10)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    mtp = mt[sel,4]/mt[sel,4][0]
    plot(mt[sel,2],(mtp-1.25)*(mtp>1.25)*2-(mtp-1.1)*(mtp>1.1),label=R2/1e6)
mtt = mt[sel,:]
sel = array([12,60])
mtp = mtp[sel]
mtt = mtt[sel,2]
plot(mtt,(mtp-1.25)*(mtp>1.25)*2-(mtp-1.1)*(mtp>1.1),'C1:')
ylabel("$\Delta$ w$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(round(R))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
ylim(-0.2,0.45)
xlim(-100,100)
savefig("1T-STDP_GABAposp5.png")



for i in range(10): close()




it = 2
xc = log10(array(Rnecks[1:-2])); xc = xc-min(xc); xc = xc/max(xc)

mes = mesR2sGprep5
R = Rdends[it]
#figure(it)
sel1 = abs(mes[:,0]-R)<1e-3
mt = mes[sel1,:]

figure(5)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    plot(mt[sel,2],mt[sel,4]/mt[sel,4][0],label=R2/1e6)
    #plot(mt[sel,2],-mt[sel,3]*(1-(mt[sel,6]>0))/max(mt[sel,3])*0.3+mt[sel,6]*0.6,'--',label=R2/1e6,c=cmapg(xc[j]))
    #plot(mes[:,1],-mes[:,2]*(1-(mes[:,4]>0))/max(mes[:,2])+mes[:,4]*0.5)
    #plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',label=R2/1e6,c=cmapg(xc[j]))

ylabel("$\Delta$ F/F$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(round(R))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
#ylim(-0.3,0.4)
savefig("1T-STDP_CA_GABAprep5.png")

#mes = mesR2sGprep5

figure(10)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    mtp = mt[sel,4]/mt[sel,4][0]
    plot(mt[sel,2],(mtp-1.25)*(mtp>1.25)*2-(mtp-1.1)*(mtp>1.1),label=R2/1e6)

ylabel("$\Delta$ w$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(round(R))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
ylim(-0.2,0.45)
xlim(-100,100)
savefig("1T-STDP_GABAprep5.png")




for i in range(10): close()




it = 2
xc = log10(array(Rnecks[1:-2])); xc = xc-min(xc); xc = xc/max(xc)

mes = mesR2sGprem5
R = Rdends[it]
#figure(it)
sel1 = abs(mes[:,0]-R)<1e-3
mt = mes[sel1,:]

figure(5)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    plot(mt[sel,2],mt[sel,4]/mt[sel,4][0],label=R2/1e6)
    #plot(mt[sel,2],-mt[sel,3]*(1-(mt[sel,6]>0))/max(mt[sel,3])*0.3+mt[sel,6]*0.6,'--',label=R2/1e6,c=cmapg(xc[j]))
    #plot(mes[:,1],-mes[:,2]*(1-(mes[:,4]>0))/max(mes[:,2])+mes[:,4]*0.5)
    #plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',label=R2/1e6,c=cmapg(xc[j]))

ylabel("$\Delta$ F/F$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(round(R))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
#ylim(-0.3,0.4)
savefig("1T-STDP_CA_GABAprem5.png")

#mes = mesR2sGprep5

figure(10)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    mtp = mt[sel,4]/mt[sel,4][0]
    plot(mt[sel,2],(mtp-1.25)*(mtp>1.25)*2-(mtp-1.1)*(mtp>1.1),label=R2/1e6)

ylabel("$\Delta$ w$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(round(R))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
ylim(-0.2,0.45)
xlim(-100,100)
savefig("1T-STDP_GABAprem5.png")



for i in range(10): close()


it = 2
xc = log10(array(Rnecks[1:-2])); xc = xc-min(xc); xc = xc/max(xc)

mes = mesRR
R = Rdends[it]
#figure(it)
sel1 = abs(mes[:,0]-R)<1e-3
mt = mes[sel1,:]

figure(5)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    plot(mt[sel,2],mt[sel,4]/mt[sel,4][0],label=R2/1e6)
    #plot(mt[sel,2],-mt[sel,3]*(1-(mt[sel,6]>0))/max(mt[sel,3])*0.3+mt[sel,6]*0.6,'--',label=R2/1e6,c=cmapg(xc[j]))
    #plot(mes[:,1],-mes[:,2]*(1-(mes[:,4]>0))/max(mes[:,2])+mes[:,4]*0.5)
    #plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',label=R2/1e6,c=cmapg(xc[j]))

ylabel("$\Delta$ F/F$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(round(R))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
#ylim(-0.3,0.4)
savefig("1T-STDP_CA_noGABA.png")

#mes = mesR2sGprep5

figure(10)
for j,R2 in enumerate(Rnecks[1:-2]):
    sel = abs(mt[:,1]-R2)<1e-3
    mtp = mt[sel,4]/mt[sel,4][0]
    plot(mt[sel,2],(mtp-1.25)*(mtp>1.25)*2-(mtp-1.1)*(mtp>1.1),label=R2/1e6)
    mtt = mt[sel,:]
    sel = array([11,60])
    if j==0: sel[1]=56
    mtp = mtp[sel]
    mtt = mtt[sel,2]
    plot(mtt,(mtp-1.25)*(mtp>1.25)*2-(mtp-1.1)*(mtp>1.1),'C'+str(j)+':')

ylabel("$\Delta$ w$")
xlabel("$\Delta t")
title("Z$_{dend}$ = "+str(round(R))+" M $\Omega$")
legend(title="R$_{neck}$ (M$\Omega$)")
ylim(-0.2,0.45)
xlim(-100,100)
savefig("1T-STDP_noGABA.png")
