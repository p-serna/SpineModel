from numpy import *
from numpy.random import randn,randint
from matplotlib.pylab import subplots
import pickle

datasetfile = "data_Fullset/Fullset.pkl"
with open(datasetfile,"rb") as f:
    data = pickle.load(f)

with open("data_Fullset/shPSD_Morphometry.pkl","rb") as f:
    datash = pickle.load(f)

def dataset(data,noise = 0.1):
    size = data[list(data.keys())[0]].shape[0]
    ies = arange(size)
    bt = {}
    for key in data.keys():
        bt[key] = data[key][ies]*clip(1.0+randn(size)*noise,0,None)

    Dmax = bt["maxDhead"]/1e3

    bt["L"] = Dmax
    bt["D"] = sqrt(4*bt["Vh"]/Dmax/pi)
    bt['AhA0'] = bt["Ah"]/(bt["D"]*bt["L"]*pi)
    
    return(bt)

def dataset_cd(data = None,cd = None,noise = 0.1,):
    if data is None:
        with open(datasetfile,"rb") as f:
            data = pickle.load(f)
    data2 = {}
    if cd is None:
        data2 = data
    if cd == 'DiS':
        sel = data['nPSD']==2
        for key in data.keys():
            data2[key] = data[key][sel]
    if cd == 'SA':
        sel = data['SA']==1
        for key in data.keys():
            data2[key] = data[key][sel]
    if cd == 'noSA':
        sel = data['SA']!= 1
        for key in data.keys():
            data2[key] = data[key][sel]
    if cd == 'SiS':
        sel = data['nPSD']==1
        for key in data.keys():
            data2[key] = data[key][sel]
    if cd == 'Sp':
        sel = (data['nPSD']>=1)
        for key in data.keys():
            data2[key] = data[key][sel]
    return(dataset(data2,noise))

def btShInh(size = 100,noise = 0.1):
    ies = randint(datash.shape[0],size=size)
    bt = datash['A'][ies]*clip(1.0+randn(size)*noise,0,None)
    
    return(bt)
def btShInhwpos(size = 100,noise = 0.1):
    ies = randint(datash.shape[0],size=size)
    bt = datash['A'][ies]*clip(1.0+randn(size)*noise,0,None)
    pos = datash['Dps'][ies]*clip(1.0+randn(size)*noise,0,None)
   
    return(bt,pos)
    
def ShInhwpos(noise = 0.1):
    size = len(datash['A'])
    bt = datash['A']*clip(1.0+randn(size)*noise,0,None)
    pos = datash['Dps']*clip(1.0+randn(size)*noise,0,None)
   
    return(bt,pos)

def btset(data,size = 100,noise = 0.1):
    ies = randint(data[list(data.keys())[0]].shape[0],size=size)
    bt = {}
    for key in data.keys():
        bt[key] = data[key][ies]*clip(1.0+randn(size)*noise,0,None)

    Dmax = bt["maxDhead"]/1e3

    bt["L"] = Dmax
    bt["D"] = sqrt(4*bt["Vh"]/Dmax/pi)
    bt['AhA0'] = bt["Ah"]/(bt["D"]*bt["L"]*pi)
    
    return(bt)

def btset_cd(data = None,cd = None,size = 100,noise = 0.1,):
    if data is None:
        with open(datasetfile,"rb") as f:
            data = pickle.load(f)
    data2 = {}
    if cd is None:
        data2 = data
    if cd == 'DiS':
        sel = data['nPSD']==2
        for key in data.keys():
            data2[key] = data[key][sel]
    if cd == 'SA':
        sel = data['SA']==1
        for key in data.keys():
            data2[key] = data[key][sel]
    if cd == 'noSA':
        sel = data['SA']!= 1
        for key in data.keys():
            data2[key] = data[key][sel]
    if cd == 'SiS':
        sel = data['nPSD']==1
        for key in data.keys():
            data2[key] = data[key][sel]
    if cd == 'Sp':
        sel = (data['nPSD']>=1)
        for key in data.keys():
            data2[key] = data[key][sel]
    return(btset(data2,size,noise))
    
    
def generateplots(spines,corplot = False, ptitle = ""):
    spkeys = list(spines.keys())
    nkeys = len(spkeys)
    nr = nkeys//3+1
    fig1, ax = subplots(nrows=nr,ncols=3,figsize=(20,6*nr))
    fig1.suptitle(ptitle,fontsize = 20)

    for i,key in enumerate(spines.keys()):
        ax[i//3,i%3].hist(spines[key],31)
        ax[i//3,i%3].set_xlabel(key)
    
    if corplot == True:
        nt = sum(arange(nkeys))
        nr = nt//3+1
        fig2, ax = subplots(nrows=nr,ncols=3,figsize=(20,6*nr))
        fig2.suptitle("Crosscor",fontsize = 20)
        k = 0
        for i in range(nkeys):
            key = spkeys[i]
            for j in range(i+1,nkeys):
                key2 = spkeys[j]
                ax[k//3,k%3].plot(spines[key],spines[key2],'.',alpha=0.7)
                ax[k//3,k%3].set_xlabel(key)
                ax[k//3,k%3].set_ylabel(key2)
                k += 1
    else:
        fig2 = ""

    
    figs = (fig1,fig2)
    return(figs)

def simulateSet(spn,tG = 500,ton = 50,toffset = 50,t_stop = 250, EL = -65,btsr = None,VDCC = array([0.,0,0])):
    if btsr is None:
        btsr =ones(7)==1
        
    dend.Ra = 250
    neck.Ra = 250
    model.E_PAS = -65
    model.soma.e_pas = model.E_PAS
    for dendp in model.dend:
        dendp.e_pas = model.E_PAS
    for sp in model.spne:
        sp.e_pas = model.E_PAS

    dendsh = model.dend[-2]
    dendc = model.dend[1]
    dendN = model.dend[-1]
    
    data = column_stack((spn["A1"],spn["A1"],spn["A2"],spn["Rn"],spn["Dss"],dis["L"],dis["D"]))
    for i in range(7):
        if btsr[i]:
            data[:,i] = data[:,i].mean()
            
    mes = zeros((nsp,9))
    vavg = zeros((int(t_stop/lb.h.dt+1),7))
    vtracs = zeros((int(t_stop/lb.h.dt+1),500))
    Ctracs = zeros((int(t_stop/lb.h.dt+1),500))
    vtracsD = zeros((int(t_stop/lb.h.dt+1),500))
    vtracsS = zeros((int(t_stop/lb.h.dt+1),500))
    for i in arange(nsp):
        NC.weight[0]  = data[i,0] *gtrA#/2
        NCN.weight[0] = data[i,1] *gtrN#*0#*0
        NCG.weight[0] = data[i,2] *gtrG#*0
        
        neck.L = 1.5
        neck.diam = diam0*sqrt(Rneck0/data[i,3])
        
        posD = data[i,4]
        
        dendsh.L = 1.0
        dendc.L = posD-0.5
        dendN.L = dendsizeL-posD-0.5
        dendsh.diam = 0.5
        dendc.diam = 0.5
        dendN.diam = 0.5
        dendN.cm = 3.5
        dendc.cm = 3.5
        
        
        sp = model.spne[0]
        # A = pi*D**2
        sp.L = data[i,5]
        sp.diam = data[i,6]
        spvol = sp(0.5).volume()

        #spineArea =  spn["Ah"][i] # um^2
        #sp.gbar_itL = 5*CaTcond/spineArea*1e-10 # S/ cm^2   not pS/um^2 as in mod file
        spineArea =  sp(0.5).area()#sp.L*sp.diam+1.8*sp.diam**2/4 # um^2
        CaTcond = 2e6# pS
        sp.gbar_itL = VDCC[0]*0.15*CaTcond/spineArea*1e-10
        sp.gbar_ca = VDCC[1]*100*CaTcond/spineArea*1e-10
        sp.gbar_it = VDCC[2]*0.15*CaTcond/spineArea*1e-10
        spineArea =  sp(0.5).area()
        
        NC.delay = toffset+ton-50
        NCN.delay = toffset+ton-50
        NCG.delay = toffset+tG#-50

        
        # ~ print(NC.weight[0])
        lb.h.finitialize(model.E_PAS)
        # ~ print(NC.weight[0])
        

        
        lb.neuron.run(t_stop)
        
        #plot(trec,vspneRec[0])    
        
        current = abs((array(vDendRec[-2])-array(vDendRec[0]))/Rdend)
        
        vtracs[:,i] = array(vspneRec[0]) 
        vtracsD[:,i] = array(vDendRec[1]) 
        vtracsS[:,i] = array(vrec) 

        vavg[:,0] += array(vspneRec[0]) 
        vavg[:,1] += array(vspneRec[0])**2
        vavg[:,2] += array(vDendRec[0]) 
        vavg[:,3] += array(vDendRec[0])**2
        vavg[:,4] += array(vrec) 
        vavg[:,5] += array(vrec)**2
        vavg[:,6] += 1
            
        cat = array(caDendRec[-1])/1e-3
        Ctracs[:,i] = cat-cat[0] 
        aG = abs(array(currentGABA)).argmax()
        aA = abs(array(currentAMPA)).argmax()
        
        mes[i,:] = [spn["Rn"][i],max(vspneRec[0]),max(vDendRec[1]),max(vrec),max(cat)-cat[0],array(currentGABA)[aG],array(currentAMPA)[aA],spvol,max(current)]
        
        #plot(trec,array(caDendRec[-1])/1e-3)
        #ylabel("[Ca] (uM)")
        #figure()
        #plot(trec,vspneRec[0])
        #break
        
    vavg[:,:5] = vavg[:,:5]/vavg[0,6]
    vavg[:,1] = sqrt(vavg[:,1]-vavg[:,0]**2)#/sqrt(vavg[0,6])
    vavg[:,3] = sqrt(vavg[:,3]-vavg[:,2]**2)#/sqrt(vavg[0,6])
    vavg[:,5] = sqrt(vavg[:,5]-vavg[:,4]**2)#/sqrt(vavg[0,6])
    return(vavg,mes,vtracs,vtracsD,vtracsS)


def ap(t,ton= 20, vb = -65,vth = -40, vp = 30, ttopeak = 1.0,tdec = 2.5,ttmax=0.1):
    ton2 = ton+ttopeak
    ton3 = ton+ttopeak+ttmax
    vtoth = (vb+(vth-vb)*(array(t)-ton)/ttopeak)*(t>ton)*(t<ton2)
    vtopeak = +(vth+(array(t)-ton2)*(vp-vth)/ttmax)*(t>ton2)*(t<ton3)
    vdec = ((vp-vb)*exp(-(array(t)-ton3)/tdec*2)+vb)*(t>ton3)

    vap = vb*(t<=ton) + vtoth + vtopeak +vdec
    return(vap)
    


def simulateSetwfAP(spn,y,t,tG = 500,ton = 50,toffset = 50,t_stop = 250, EL = -65,btsr = None,VDCC = array([0.,0,0]),tAP = 100):
    if btsr is None:
        btsr =ones(7)==1
        
    neck.Ra = 250
    model.E_PAS = -65
    model.soma.e_pas = model.E_PAS
    for dendp in model.dend:
        dendp.e_pas = model.E_PAS
    for sp in model.spne:
        sp.e_pas = model.E_PAS

    dendsh = model.dend[-2]
    dendc = model.dend[1]
    dendN = model.dend[-1]
    
    data = column_stack((spn["A1"],spn["A1"],spn["A2"],spn["Rn"],spn["Dss"],spn["L"],spn["D"]))
    for i in range(7):
        if ~btsr[i]:
            data[:,i] = data[:,i].mean()
    
    print(data.std(axis=0))
            
    y.play_remove()
    y.play(dendsh(0.5)._ref_v,t,True)
    y.play(dendc(0.5)._ref_v,t,True)
    y.play(model.soma(0.5)._ref_v,t,True)

    mes = zeros((nsp,9))
    vavg = zeros((int(t_stop/lb.h.dt+1),7))
    vtracs = zeros((int(t_stop/lb.h.dt+1),500))
    Ctracs = zeros((int(t_stop/lb.h.dt+1),500))
    vtracsD = zeros((int(t_stop/lb.h.dt+1),500))
    vtracsS = zeros((int(t_stop/lb.h.dt+1),500))
    for i in arange(nsp):
        NC.weight[0]  = data[i,0] *gtrA#/2
        NCN.weight[0] = data[i,1] *gtrN#*0#*0
        NCG.weight[0] = data[i,2] *gtrG#*0
        
        neck.L = 1.5
        neck.diam = diam0*sqrt(Rneck0/data[i,3])
        
        posD = data[i,4]
        
        dendsh.L = 1.0
        dendc.L = posD-0.5
        dendN.L = dendsizeL-posD-0.5
        dendsh.diam = 0.5
        dendc.diam = 0.5
        dendN.diam = 0.5
        dendN.cm = 3.5
        dendc.cm = 3.5
        
        
        sp = model.spne[0]
        # A = pi*D**2
        sp.L = data[i,5]
        sp.diam = data[i,6]
        spvol = sp(0.5).volume()

        #spineArea =  spn["Ah"][i] # um^2
        #sp.gbar_itL = 5*CaTcond/spineArea*1e-10 # S/ cm^2   not pS/um^2 as in mod file
        spineArea =  sp(0.5).area()#sp.L*sp.diam+1.8*sp.diam**2/4 # um^2
        CaTcond = 2e6# pS
        sp.gbar_itL = VDCC[0]*0.15*CaTcond/spineArea*1e-10
        sp.gbar_ca = VDCC[1]*100*CaTcond/spineArea*1e-10
        sp.gbar_it = VDCC[2]*0.15*CaTcond/spineArea*1e-10
        spineArea =  sp(0.5).area()
        
        NC.delay = toffset+ton-50
        NCN.delay = toffset+ton-50
        NCG.delay = toffset+tG#-50

        ta = concatenate(([0],linspace(tAP-10,tAP+40,10000),[t_stop]))

        y = y.from_python(ap(ta,tAP))
        t = t.from_python(ta)
        
        # ~ print(NC.weight[0])
        lb.h.finitialize(model.E_PAS)
        # ~ print(NC.weight[0])
        

        
        lb.neuron.run(t_stop)
        
        #plot(trec,vspneRec[0])    
        
        current = abs((array(vDendRec[-2])-array(vDendRec[0]))/Rdend)
        
        vtracs[:,i] = array(vspneRec[0]) 
        vtracsD[:,i] = array(vDendRec[1]) 
        vtracsS[:,i] = array(vrec) 

        vavg[:,0] += array(vspneRec[0]) 
        vavg[:,1] += array(vspneRec[0])**2
        vavg[:,2] += array(vDendRec[0]) 
        vavg[:,3] += array(vDendRec[0])**2
        vavg[:,4] += array(vrec) 
        vavg[:,5] += array(vrec)**2
        vavg[:,6] += 1
            
        cat = array(caDendRec[-1])/1e-3
        Ctracs[:,i] = cat-cat[0] 
        aG = abs(array(currentGABA)).argmax()
        aA = abs(array(currentAMPA)).argmax()
        
        mes[i,:] = [spn["Rn"][i],max(vspneRec[0]),max(vDendRec[1]),max(vrec),max(cat)-cat[0],array(currentGABA)[aG],array(currentAMPA)[aA],spvol,max(current)]
        
        #plot(trec,array(caDendRec[-1])/1e-3)
        #ylabel("[Ca] (uM)")
        #figure()
        #plot(trec,vspneRec[0])
        #break
        
    vavg[:,:5] = vavg[:,:5]/vavg[0,6]
    vavg[:,1] = sqrt(vavg[:,1]-vavg[:,0]**2)#/sqrt(vavg[0,6])
    vavg[:,3] = sqrt(vavg[:,3]-vavg[:,2]**2)#/sqrt(vavg[0,6])
    vavg[:,5] = sqrt(vavg[:,5]-vavg[:,4]**2)#/sqrt(vavg[0,6])
    return(vavg,mes,vtracs,vtracsD,vtracsS,Ctracs)
