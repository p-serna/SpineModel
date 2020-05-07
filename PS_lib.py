from neuron import *
load_mechanisms('mod/')
h('objref nil')

def parameters(model):
    # Passive properties
    model.CM = 1.0
    model.RM = 7000.0
    model.RA = 150.0 
    model.E_PAS = -70
    model.temperature = 37

    # Active properties
    model.Ek = -90
    model.Ena = 60
    model.Eca = 140
    
    model.gna_axon = 1000
    model.gkv_axon = 100
    
    model.gna_soma = 1000
    model.gkv_soma = 100 
    model.gkm_soma = 2.2 
    model.gkca_soma = 3 
    model.gca_soma = 0.5 
    model.git_soma = 0.0003 
    
    model.gna_dend = 80
    #model.gna_dend_hotSpot = 600
    model.gkv_dend = 3
    model.gkm_dend = 1
    model.gkca_dend = 3
    model.gca_dend = 0.5
    model.git_dend = 0.00015 
    model.gh_dend = 0
    
    model.NMDAlist = []
    model.ncNMDAlist = []
    model.GABAlist = []
    model.ncGABAlist = []
    model.stims = []
    model.dt = 0.05
    

class loadNeuron(object):
    def __init__(self,hocfile="Basic.hoc",axon=False):
        h('xopen("'+hocfile+'")')
        parameters(self)
        self._topol()
        if axon:
            self._axon()
        self._biophys()
        #self.topology = h.topology()

    def _topol(self):            
        self.soma = h.soma
        self.dend = []
        self.comp = {}
        for sec in h.allsec():
            self.dend.append(sec)
            self.comp[sec.hname()] = sec
        # Dropping Soma 
        self.dend.pop() 
        self.spne, self.neck = [], []
        
        
    def _axon(self):
        self.axon = h.Section(name='axon')
        self.axon.L = 300
        self.axon.diam = 1
        self.axon.connect(self.soma,1,0)
        self.comp['axon'] = self.axon 
    
    def _biophys(self):
        for sec in h.allsec():
            sec.cm = self.CM
            sec.insert('pas')
            sec.e_pas = self.E_PAS
            sec.g_pas = 1.0/self.RM
            sec.Ra = self.RA
            
    def addDend(self,name=None,locus=-1,L=2.0,D=1.0,ilocus=0.5):
        if name is None:
            name = "Dend"+str(len(self.dend)).zfill(3)
        newdend = h.Section(name=name)
        if isinstance(locus,str):
            #pass
            for sec in h.allsec():
               if locus==sec.name(): 
                   newdend.connect(sec(ilocus))
                   break
        else:
            if locus == -1:
                newdend.connect(self.soma(ilocus))
            elif locus<len(self.dend):
                dendP = self.dend[locus]
                newdend.connect(dendP(ilocus))
        newdend.cm = self.CM
        newdend.insert('pas')
        newdend.e_pas = self.E_PAS
        newdend.g_pas = 1.0/self.RM
        newdend.Ra = self.RA
        newdend.L = L
        newdend.diam = D
        
        self.dend.append(newdend)
        self.comp[name] = newdend 

    def addSpne(self,name=None,locus=-1,ilocus=0.5,L=1.0,D=1.0,Lneck=1.2,Dneck=0.1):
        if name is None:
            name = "Spne"+str(len(self.spne)).zfill(3)
        newspne = h.Section(name=name)
        newneck = h.Section(name="N"+name)
        if isinstance(locus,str):
            #pass
            for sec in h.allsec():
               if locus==sec.name(): 
                   newneck.connect(sec(ilocus))
                   newspne.connect(newneck(1))
                   break
        else:
            if locus == -1:
                newneck.connect(self.soma(ilocus))
                newspne.connect(newneck(1))
            elif locus<len(self.dend):
                dendP = self.dend[locus]
                newneck.connect(dendP(ilocus))
                newspne.connect(newneck(1))
        newspne.cm = self.CM
        newspne.insert('pas')
        newspne.e_pas = self.E_PAS
        newspne.g_pas = 1.0/self.RM
        newspne.Ra = self.RA
        newspne.L = L
        newspne.diam = D
        newneck.cm = self.CM
        newneck.insert('pas')
        newneck.e_pas = self.E_PAS
        newneck.g_pas = 1.0/self.RM
        newneck.Ra = self.RA
        newneck.L = Lneck
        newneck.diam = Dneck
                        
        self.spne.append(newspne)
        self.neck.append(newneck)
        self.comp[name] = newspne 
        self.comp[name+'_neck'] = newneck 







def init_active(model, axon=False, soma=False, dend=True, dendNa=False,
                dendCa=False,spne=False):
    if axon:
        model.axon.insert('na'); model.axon.gbar_na = model.gna_axon
        model.axon.insert('kv'); model.axon.gbar_kv = model.gkv_axon
        model.axon.ena = model.Ena
        model.axon.ek = model.Ek

    if soma:
        model.soma.insert('na'); model.soma.gbar_na = model.gna_soma
        model.soma.insert('kv'); model.soma.gbar_kv = model.gkv_soma
        model.soma.insert('km'); model.soma.gbar_km = model.gkm_soma
        model.soma.insert('kca'); model.soma.gbar_kca = model.gkca_soma
        model.soma.insert('ca'); model.soma.gbar_ca = model.gca_soma
        model.soma.insert('it'); model.soma.gbar_it = model.git_soma
        #model.soma.insert('cad');
        model.soma.ena = model.Ena
        model.soma.ek = model.Ek
        model.soma.eca = model.Eca

    if dend:
        for d in model.dend:
            d.insert('na'); d.gbar_na = model.gna_dend*dendNa
            d.insert('kv'); d.gbar_kv = model.gkv_dend
            d.insert('km'); d.gbar_km = model.gkm_dend
            d.insert('kca'); d.gbar_kca = model.gkca_dend
            d.insert('ca'); d.gbar_ca = model.gca_dend*dendCa
            d.insert('it'); d.gbar_it = model.git_dend*dendCa
            #d.insert('cad')
            d.ena = model.Ena
            d.ek = model.Ek
            d.eca = model.Eca
    if spne:
        for s in model.spne:
            #s.insert('na'); s.gbar_na = model.gna_dend*dendNa
            #s.insert('kv'); s.gbar_kv = model.gkv_dend
            #s.insert('km'); s.gbar_km = model.gkm_dend
            #s.insert('kca'); s.gbar_kca = 0*model.gkca_dend
            s.insert('cal_ion')
            s.insert('caqPS'); s.pcaqbar_caqPS = 6.0e-6      
            # N - type
            #s.insert('can'); s.pbar_can = 1.0e-5      
            s.insert('canPS'); s.pbar_canPS = 1.0e-5      
            
            # (HVA) L - type
            #s.insert('caL'); s.pbar_caL = 6.7e-6      
            s.insert('caLPS'); s.pbar_caLPS = 6.7e-6     
            
            # (Cav1.3) L - type
            #s.insert('caL13'); s.pcaLbar_caL13 = 1.7e-6
            s.insert('caL13PS'); s.pbar_caL13PS = 1.7e-6     
            s.insert('cad'); 
            #s.ca = 0.0023
            #s.insert('cad')
            #s.ena = model.Ena
            #s.ek = model.Ek
            #s.eca = model.Eca



def add_somaStim(model, p=0.5, onset=20, dur=1, amp=0):
    model.stim = h.IClamp(model.soma(p))
    model.stim.delay = onset
    model.stim.dur = dur
    model.stim.amp = amp    # nA
    

def add_NMDAsyns(model, locs=[[0, 0.5]], gmax=0.5, tau1=2, tau2=20):
    gmax = gmax/1000.   # Set in nS and convert to muS
    for loc in locs:
        NMDA = h.Exp2SynNMDA(float(loc[1]), sec=model.spne[int(loc[0])]) 
        NMDA.tau1 = tau1
        NMDA.tau2 = tau2
        NC = h.NetCon(h.nil, NMDA, 0, 0, gmax)
        x = float(loc[1])
        model.NMDAlist.append(NMDA)
        model.ncNMDAlist.append(NC)   

def add_GABAsyns(model, spne=False, locs=[[0, 0.5]], gmax=0.5, tau1=0.1, tau2=4,
                     rev=-80):
    gmax = gmax/1000.   # Set in nS and convert to muS
    if spne:
        for loc in locs:
            GABA = h.Exp2Syn(float(loc[1]), sec=model.spne[int(loc[0])]) 
            GABA.tau1 = tau1
            GABA.tau2 = tau2
            GABA.e = rev
            stimI = h.NetStim(); stimI.number = 1
            NC = h.NetCon(stimI, GABA, 0, 0, gmax)
            model.GABAlist.append(GABA)
            model.ncGABAlist.append(NC)
            model.stims.append(stimI)
            
            #model.ncstimlist.append(NC)
    else:
        for loc in locs:
            GABA = h.Exp2Syn(float(loc[1]), sec=model.dend[int(loc[0])]) 
            GABA.tau1 = tau1
            GABA.tau2 = tau2
            GABA.e = rev
            NC = h.NetCon(h.nil, GABA, 0, 0, gmax)
            model.GABAlist.append(GABA)
            model.ncGABAlist.append(NC)
    return(GABA,NC)
    

def add_GABAsynscomp(model,comp, loc=0.5, gmax=0.5, tau1=0.1, tau2=4,
                     rev=-70):
    gmax = gmax/1000.   # Set in nS and convert to muS

    GABA = h.Exp2Syn(float(loc), sec=comp) 
    GABA.tau1 = tau1
    GABA.tau2 = tau2
    GABA.e = rev
    stimI = h.NetStim(); stimI.number = 1
    NC = h.NetCon(stimI, GABA, 0, 0, gmax)    
    model.GABAlist.append(GABA)
    model.ncGABAlist.append(NC)
    model.stims.append(stimI)
    return(GABA,NC)
    
