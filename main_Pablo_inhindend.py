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
lb.add_GABAsyns(model, locs=[[1,0.3]], gmax=gmaxG,tau1=1.5,tau2=20.0)  

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
sp.gbar_itL = 5*CaTcond/spineArea*1e-10 # S/ cm^2   not pS/um^2 as in mod file

NMDA.tau1 = 2.0
NMDA.tau2 = 30.0
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
    

tdels = 10**linspace(-0.5,log10(50),30)
tdels2 = hstack((-tdels[arange(len(tdels)-1,-1,-1)],0,tdels))
t_stop = 2000
ton = 1100

NCG.weight[0] = 30*35/spineArea*1e-10*1e4
sp.gbar_itL = 5*CaTcond/spineArea*1e-10 # S/ cm^2   not pS/um^2 as in mod file
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
ylabel("Ca$^{2+}$ [ ($\Delta$ F)/F ]")
#legend()
ylim(0.82,1.02)
twinx()
plot(mes[:,0],(mes[:,1])/max(mes[:,1]),"C1--",label="Voltage ")
ylim(0.70,1.05)
ylabel("Voltage [$\Delta$ V/V]")
sel = mes[:,1]/max(mes[:,1])==1
tp = mes[sel,0]
sel = tp>0 
vlines(tp[sel][0],0,1.2,'k',linestyles="dashed")

savefig("0T-DEND-Inhvsexc.png")

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
tdels2 = hstack((-tdels[arange(len(tdels)-1,-1,-1)],0,tdels))

#-------------- VArying Neck resistance

mes = []
mp = []
Rnecks = array([1,10,20,50,100,200,500,1000,10000])*1e6

t_stop = 2000
ton = 1100
lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)

NCG.weight[0] = 30*35/spineArea*1e-10*1e4
sp.gbar_itL = 5*CaTcond/spineArea*1e-10 # S/ cm^2   not pS/um^2 as in mod file
NCN.weight[0] = 50*0.0002/1e3

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
savetxt("DwRneckRdend36MOhm.dat",mes)

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
    plot(mt[sel,1],mt[sel,3]/max(mt[sel,3]),label=R2/1e6,c=cmapf(xc[j]))
        
xlabel("Delta t")
ylabel("Calcium [ $\Delta$ F/F ]")
    #legend()    
legend(title="R$_{neck}$ (M$\Omega$)")

twinx()
for j,R2 in enumerate(Rnecks[1:-1]):
    sel = abs(mt[:,0]-R2)<1e-3
    plot(mt[sel,1],mt[sel,2]/max(mt[sel,2]),'--',label=R2/1e6,c=cmapg(xc[j]))

ylabel("$\Delta$ V/V$_{max}$")

title("Z$_{dend}$ = "+str(int(Rdend/1e6))+" M $\Omega$")

savefig("zdend"+str(int(Rdend/1e6)).zfill(4)+"mOhm.png")


# Varying both neck and dendritic resistance.

mes = []
mp = []
#Rnecks = 10**linspace(log10(1e6),log10(1e9),21)

Rnecks = array([1,10,20,50,100,200,500,1000,2000,10000])*1e6
Rdends0 = array([70,100,200,300,400,500,750,1000,1500,2000])*1e6

Rnecks = array([50,500,1000])*1e6
Rdends0 = array([70,200,500,1000])*1e6



Rs = [[i,j] for i in Rnecks for j in Rdends0]

t_stop = 2000
ton = 1100
lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)


NCG.weight[0] = 30*35/spineArea*1e-10*1e4
sp.gbar_itL = 20*CaTcond/spineArea*1e-10 # S/ cm^2   not pS/um^2 as in mod file
NCN.weight[0] = 50*0.0002/1e3


dendt = model.dend[1]
dendt.L = 50.0
Rmin = 70.84852670568594e6 #177.12856043961426


mes = []
mp = []
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
    print(Rdend2)
    NC.delay = 3050
    NCN.delay = 3050
    NCG.delay = 3050.0
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    mCat =  caDendRec[3][-1]
    mV = vspneRec[0][-1]
    
    mp.append([Rdend,Rneck,mV,mCat])
    for tdel in tdels2:
        NC.delay = ton-50
        NCN.delay = ton-50
        NCG.delay = ton+tdel-50.0

        lb.h.finitialize(model.E_PAS)
        lb.neuron.run(t_stop)
        mes.append([Rdend2,Rneck,tdel,max(vspneRec[0])-mV,max(caDendRec[3])-mCat])
    print(Rdend/1e6,Rneck/1e6)


mes = array(mes)

mesR2s = mes

savetxt("DwR2As.dat",mesR2s)

cm = 'Blues'
cm2 = 'Oranges'


#--------------- Some plots: A mejorar

figure(5)
xc = log10(array(Rnecks)); xc = xc-min(xc); xc = xc/max(xc)
 xc = array([0.25,0.8])

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
    for j,R2 in enumerate(Rnecks[:-1]):
        sel = abs(mt[:,1]-R2)<1e-3
        plot(mt[sel,2],mt[sel,4]/max(mt[sel,4]),label=R2/1e6,c=cmapf(xc[j]))
            
    xlabel("Delta t")
    ylabel("Calcium [ $\Delta$ F/F ]")
    legend(title="R$_{Neck}$ (M$\Omega$)")    
        
    twinx()
    for j,R2 in enumerate(Rnecks[:-1]):
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



Rnecks = array([50,500])*1e6
Rdends0 = array([100,500])*1e6



Rs = [[i,j] for i in Rnecks for j in Rdends0]

t_stop = 2000
ton = 1100
lb.add_somaStim(model,p=0.5,onset=ton+1e4,dur=2, amp =0)


NCG.weight[0] = 30*35/spineArea*1e-10*1e4
sp.gbar_itL = 20*CaTcond/spineArea*1e-10 # S/ cm^2   not pS/um^2 as in mod file
NCN.weight[0] = 50*0.0002/1e3


dendt = model.dend[1]
dendt.L = 50.0
Rmin = 70.84852670568594e6 #177.12856043961426


mes = []
mes2 = []
mp = []
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
    print(Rdend2)
    NC.delay = 3050
    NCN.delay = 3050
    NCG.delay = 3050.0
    lb.h.finitialize(model.E_PAS)
    lb.neuron.run(t_stop)
    mCat =  caDendRec[3][-1]
    mV = vspneRec[0][-1]
    
    mp.append([Rdend,Rneck,mV,mCat])
    #for tdel in tdels2[[0,13,40,67,-1]]:
    for tdel in tdels2:
    #~ for tdel in [tdels2[57]]:
        m =[]
        for i in range(200):
            NC.delay = ton-50
            NCN.delay = ton-50
            NCG.delay = ton+tdel-50.0
            NMDA.tau2 = rand()*60+20
            nampa = random.lognormal(mean=3.9123954,sigma=0.40)
            ngaba = random.lognormal(mean=3.9123954,sigma=0.40)

            Rneckt = Rneck
            neck.Ra = Rneckt/neck.L*(neck.diam*1e-6/2.0)**2*pi*100/1e-6
            NCG.weight[0] = 30*ngaba/spineArea*1e-10*1e4
            #NCN.weight[0] = 50*0.0002/1e3
            NC.weight[0]  =  15*nampa/1e6

            lb.h.finitialize(model.E_PAS)
            lb.neuron.run(t_stop)
            #plot(trec,caDendRec[3],'C0-',alpha=0.4)
            m.append([max(vspneRec[0])-mV,max(caDendRec[3])-mCat])
            if i%10==0: print(i,tdel)
        m = array(m)
        h = histogram(m[:,1],bins=30)
        sel1 = arange(len(h[0]))[cumsum(h[0])/sum(h[0])>0.25][0]
        sel2 = arange(len(h[0]))[cumsum(h[0])/sum(h[0])>0.75][0]
        dm1 = h[1][[sel1,sel2]]
        h = histogram(m[:,0],bins=30)
        sel1 = arange(len(h[0]))[cumsum(h[0])/sum(h[0])>0.25][0]
        sel2 = arange(len(h[0]))[cumsum(h[0])/sum(h[0])>0.75][0]
        dm0 = h[1][[sel1,sel2]]
        m0 = mean(m,axis=0)
        mes.append([Rdend2,Rneck,tdel,m0[0],m0[1],dm0[0],dm0[1],dm1[0],dm1[1]])
        m =[]
        for i in range(200):
            NC.delay = ton-50
            NCN.delay = ton-50
            NCG.delay = ton+tdel-50.0
            NMDA.tau2 = rand()*60+20
            nampa = 50 #random.lognormal(mean=3.6123954,sigma=0.10)
            ngaba = 50 #random.lognormal(mean=3.9123954,sigma=0.10)
            
            Rneckt = random.lognormal(mean=log(Rneck),sigma=0.40)
            neck.Ra = Rneckt/neck.L*(neck.diam*1e-6/2.0)**2*pi*100/1e-6

            NCG.weight[0] = 30*ngaba/spineArea*1e-10*1e4
            #NCN.weight[0] = 50*0.0002/1e3
            NC.weight[0]  =  15*nampa/1e6

            lb.h.finitialize(model.E_PAS)
            lb.neuron.run(t_stop)
            #plot(trec,caDendRec[3],'C0-',alpha=0.4)
            m.append([max(vspneRec[0])-mV,max(caDendRec[3])-mCat])
            if i%10==0: print(i,tdel)
        m = array(m)
        h = histogram(m[:,1],bins=30)
        sel1 = arange(len(h[0]))[cumsum(h[0])/sum(h[0])>0.25][0]
        sel2 = arange(len(h[0]))[cumsum(h[0])/sum(h[0])>0.75][0]
        dm1 = h[1][[sel1,sel2]]
        h = histogram(m[:,0],bins=30)
        sel1 = arange(len(h[0]))[cumsum(h[0])/sum(h[0])>0.25][0]
        sel2 = arange(len(h[0]))[cumsum(h[0])/sum(h[0])>0.75][0]
        dm0 = h[1][[sel1,sel2]]
        m0 = mean(m,axis=0)
        mes2.append([Rdend2,Rneck,tdel,m0[0],m0[1],dm0[0],dm0[1],dm1[0],dm1[1]])

    print(Rdend/1e6,Rneck/1e6)

mes = array(mes)
mes2 = array(mes2)

mesRs = mes
mesNeck = mes2

#savetxt("dinhexRs10.dat",mes)
savetxt("dinhexRs20.dat",mes)

#savetxt("dinhexneck10.dat",mes2)
savetxt("dinhexneck20.dat",mes2)

#mesonlyT = mes
#----- Este fue para el otro
mt = mes
for i in range(10): close()

fill_between(mt[:,2],mt[:,8]/max(mt[:,4]),mt[:,7]/max(mt[:,4]),color=cmapf(xc[j]),alpha=0.5)
plot(mt[:,2],mt[:,4]/max(mt[:,4]),'C0')

ylabel("$\Delta$ F/F")
xlabel("$\Delta$ t(ms)")
twinx()
fill_between(mt[:,2],mt[:,5]/max(mt[:,3]),mt[:,6]/max(mt[:,3]),color='C1',alpha=0.5)
plot(mt[:,2],mt[:,3]/max(mt[:,3]),'C1--')

ylabel("$\Delta$ V/V")
#savefig("FlucttauNMDA.png")
#savefig("FlucttauRs10.png")
savefig("FlucttauRs20.png")
#savefig("FluctNeck10.png")

mt = mes2
for i in range(10): close()

fill_between(mt[:,2],mt[:,8]/max(mt[:,4]),mt[:,7]/max(mt[:,4]),color=cmapf(xc[j]),alpha=0.5)
plot(mt[:,2],mt[:,4]/max(mt[:,4]),'C0')

ylabel("$\Delta$ F/F")
xlabel("$\Delta$ t(ms)")
twinx()
fill_between(mt[:,2],mt[:,5]/max(mt[:,3]),mt[:,6]/max(mt[:,3]),color='C1',alpha=0.5)
plot(mt[:,2],mt[:,3]/max(mt[:,3]),'C1--')

ylabel("$\Delta$ V/V")
#savefig("FlucttauNMDA.png")
#savefig("FlucttauRs10.png")
#savefig("FluctNeck10.png")
savefig("FluctNeck20.png")

