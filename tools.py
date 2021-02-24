# set of tools that are used everywhere
import numpy as np
from scipy.integrate import simps

def showmessage(message, length = 40):
    print(message+"."*(length-len(message)))

def itertis(t1,t2,A1,A2):
    te = t1*t2/(t1-t2)
    numtt = (t2/t1)**(te/t1)-(t2/t1)**(te/t2)
    t1 = t2 + numtt*A1
    t2 = A2/np.log(t1/t2)*(t1-t2)/t1
    return(t1,t2)
    
def gettimes(As,t1 = 10,t2 = 1, n = 1000,tol = 1e-9):
    A1, A2 = As
    t1a, t2a = t1, t2
    for i in range(n):
        t1,t2 = itertis(t1,t2,A1,A2)
    
        if (t1-t1a)**2/t1**2+(t2-t2a)**2/t2**2<tol:
            #print(i)
            break
        t1a, t2a = t1, t2
        #print(t1,t2)
    if i>=n: 
        print('No convergence?')
        
    return(t1,t2)
    
def get_atimes(vtt,ton = 3000, dt = 0.05):
    ampv = vtt.max(axis=0)
    trec = np.arange(vtt.shape[0])*dt
    intv = np.array(list(map(lambda x: simps(x,trec),vtt.transpose())))
    tt = intv/ampv
    sel = vtt.argmax(axis=0)
    trise = (sel-ton)*dt
    times = np.array(list(map(gettimes,np.column_stack((tt,trise)))))
    return(times, tt, trise)

def get_FWHMtimes(vtt,v0 = None,dt =0.05):
    if v0 is None:
        v0 = vtt[0]
    ampv = vtt.max(axis=0,keepdims=True)-v0
    ampv2 = ampv/2+v0
    sel = ((vtt-ampv2)>0).transpose() 
    times = np.array(list(map(lambda x: [np.arange(len(x))[x][0],np.arange(len(x))[x][-1]],sel)))
    times = dt*(times[:,1]-times[:,0])
    return(times)
    
def get_p0p1times(vtt,p0=0.2,p1= 0.8, dt = 0.05):
    vmin = vtt.min()
    sel = vtt.argmax(axis=0)
    ampv = vtt.max(axis=0,keepdims=True)-vmin
    idcs = np.repeat(np.arange(vtt.shape[0]).reshape(vtt.shape[0],1),vtt.shape[1],1)
    sel = (vtt-vmin<ampv*p1)*(vtt-vmin>ampv*p0)*(idcs<sel)
    
    trise = (sel).sum(axis=0)*dt
    return(trise)

def get_postp0p1times(vtt,p0=0.2,p1= 0.8, dt = 0.05):
    vmin = vtt.min()
    sel = vtt.argmax(axis=0)
    ampv = vtt.max(axis=0,keepdims=True)-vmin
    idcs = np.repeat(np.arange(vtt.shape[0]).reshape(vtt.shape[0],1),vtt.shape[1],1)
    sel = (vtt-vmin<ampv*p1)*(vtt-vmin>ampv*p0)*(idcs>sel)
    
    tdecay = (sel).sum(axis=0)*dt
    return(tdecay)
    
def getint(x):
    xc = x*1.0
    xc.sort()
    sh = xc.shape[0]
    xmed = xc[sh//2]
    s0= int(sh*(1-.6827)/2)
    s1 = sh-s0-1
    x0 = xc[s0]
    x1 = xc[s1]
    s0b= int(sh*(1-.95)/2)
    s1b = sh-s0b
    try:
        x0b = xc[s0b]
        x1b = xc[s1b]
    except:
        x0b, x1b = xc[0], xc[-1]
    return((xmed,x0,x1,x0b,x1b))
    
def getintp(x,p = .95):
    xc = x*1.0
    xc.sort()
    sh = xc.shape[0]
    xmed = xc[sh//2]
    s0= int(sh*(1-p)/2)
    s1 = sh-s0-1
    x0 = xc[s0]
    x1 = xc[s1]
    return((xmed,x0,x1))
    
def getq(x,p = .95):
    xc = x*1.0
    xc.sort()
    sh = xc.shape[0]
    s0= int(sh*p)
    x0 = xc[s0]
    return(x0)

def getp(x,xc):
    xt = x*1.0
    xt.sort()    
    return((xt>xc).sum()/len(xt))
