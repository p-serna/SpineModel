

from scipy.integrate import quad
from numpy import *
from matplotlib.pylab import *
import matplotlib as mpl
from matplotlib.transforms import Bbox


import PIL.Image as Image
import pickle
import os


assert len(sys.argv)==2, 'No correct number of arguments'


condition = sys.argv[1]
folderstore = '/mnt/data/spinemodel/conditions_strongexc/'
folderoutput = '/mnt/data/spinemodel/conditions_strongexc/'+condition+'/'

if not os.path.exists(folderoutput):
    os.makedirs(folderoutput)
    



rc('text',usetex=True)
rc('text.latex',preamble='''\\usepackage{amssymb}\n\\usepackage{siunitx}\n\DeclareSIUnit\Molar{\\textsc{m}}\n
\def\\V{\\textrm{V}}\n
\def\\A{\\textrm{A}}
\def\\C{\\textrm{C}}
\def\\R{\\textrm{R}}
\def\\t{\\textrm{t}}
''')

rcParams['savefig.pad_inches'] = 0


# In[4]:


redcolor = '#dd1c77' # actually magenta


# In[5]:


SMALL_SIZE = 13
MEDIUM_SIZE = 16
BIGGER_SIZE = 20

rc('font',size = BIGGER_SIZE,family = 'Arial')
rc('axes',titlesize = BIGGER_SIZE)
rc('axes',labelsize = BIGGER_SIZE)
rc('xtick',labelsize = BIGGER_SIZE)
rc('ytick',labelsize = BIGGER_SIZE)
rc('legend',handlelength= 1.0,fontsize = MEDIUM_SIZE)
rc('figure',titlesize = BIGGER_SIZE)


# In[6]:


alphabet = []
for letter in range(65, 91):
    alphabet.append(chr(letter))
    
print(alphabet)


# In[7]:


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
    x0b = xc[s0b]
    x1b = xc[s1b]
    return((xmed,x0,x1,x0b,x1b))


# In[8]:


def numbering_panels(axs,pos = None,labels=alphabet):
    if pos is None:
        pos = zeros((len(axs),2))
        pos[:,1] = 1-pos[:,1]
        
    for i,ax in enumerate(axs):
        ax.text(pos[i,0],pos[i,1],labels[i],horizontalalignment='right',verticalalignment='top', transform=ax.transAxes)
    return


# In[9]:


def plot_trace(data,t,ax=None,c='C0',band= None,label= None,linestyle='-'):
    if ax is None:
        ax = gca()
    vtracso = data*1.0
    vtracso.sort(axis=1)
    sh = vtracso.shape[1]
    nmed = sh//2
    nl1, nv1 = int(sh*0.16),int(sh*(1-0.16))
    nl2, nv2 = int(sh*0.025),int(sh*(1-0.025))

    if label is None:
        ax.plot(t,vtracso[:,nmed],c=c,linestyle=linestyle)
    else:
        ax.plot(t,vtracso[:,nmed],c=c,linestyle=linestyle,label=label)
        
    if band == 0:
        pass
    elif band == 1:
        ax.fill_between(t,vtracso[:,nl1],vtracso[:,nv1],color=c,alpha=0.5)
    elif band ==2:
        ax.fill_between(t,vtracso[:,nl2],vtracso[:,nv2],color=c,alpha=0.5)
    else:
        ax.fill_between(t,vtracso[:,nl1],vtracso[:,nv2],color=c,alpha=0.2)
        ax.fill_between(t,vtracso[:,nl1],vtracso[:,nv1],color=c,alpha=0.5)
    return(ax)

# Scale bars
def scalebar(ax,x0,y0,dx,dy,xlab = '', ylab = '', color = 'k'):
    xs = x0+linspace(0,dx,4)
    ax.vlines(x0,y0,y0+dy,color = color)
    ax.plot(xs,xs*0+y0,c = color)
    ax.text(x0,y0+dy*.3,ylab,horizontalalignment="right")   
    ax.text(x0+dx/2,y0-dy*0.02,xlab,verticalalignment="top",horizontalalignment="center")     

def running_mean(x, N, padding = "valid"):
    if padding =="same":
        N2 = N//2
        Nf = float(N)
        cumsumt = cumsum(concatenate((zeros(1),x,zeros(N-1))))
        runmean = (cumsumt[N:] - cumsumt[:-N]) / Nf
        runmean[-N+1:] = runmean[-N+1:]*Nf/(arange(Nf-1,0,-1))
    elif padding =="valid":
        cumsumt = cumsum(insert(x, 0, 0))
        runmean = (cumsumt[N:] - cumsumt[:-N]) / float(N)
    return(runmean)
def running_std(x, N, padding = "valid"):
    if padding =="same":
        N2 = N//2
        Nf = float(N)
        cumsumt = cumsum(concatenate((zeros(1),x,zeros(N-1))))
        runmean = (cumsumt[N:] - cumsumt[:-N]) / Nf
        runmean[-N+1:] = runmean[-N+1:]*Nf/(arange(Nf-1,0,-1))
        cumsumt = cumsum(concatenate((zeros(1),x*x,zeros(N-1))))
        runstd = (cumsumt[N:] - cumsumt[:-N]) / Nf
        runstd[-N+1:] = runmean[-N+1:]*Nf/(arange(Nf-1,0,-1))
        runstd = sqrt(runstd-runmean**2)
    elif padding =="valid":
        cumsumt = cumsum(insert(x, 0, 0))
        runmean = (cumsumt[N:] - cumsumt[:-N]) / float(N)
        
        cumsumt = cumsum(insert(x**2, 0, 0))
        runstd = (cumsumt[N:] - cumsumt[:-N]) / float(N)
        runstd = sqrt(runstd-runmean**2)
    return(runstd)


# In[10]:


def gradient_image(ax, extent, direction=0.3, cmap_range=(0, 1), **kwargs):
    """
    Draw a gradient image based on a colormap.

    Parameters
    ----------
    ax : Axes
        The axes to draw on.
    extent
        The extent of the image as (xmin, xmax, ymin, ymax).
        By default, this is in Axes coordinates but may be
        changed using the *transform* kwarg.
    direction : float
        The direction of the gradient. This is a number in
        range 0 (=vertical) to 1 (=horizontal).
    cmap_range : float, float
        The fraction (cmin, cmax) of the colormap that should be
        used for the gradient, where the complete colormap is (0, 1).
    **kwargs
        Other parameters are passed on to `.Axes.imshow()`.
        In particular useful is *cmap*.
    """
    phi = direction * np.pi / 2
    v = np.array([np.cos(phi), np.sin(phi)])
    X = np.array([[v @ [1, 0], v @ [1, 1]],
                  [v @ [0, 0], v @ [0, 1]]])
    a, b = cmap_range
    X = a + (b - a) / X.max() * X
    im = ax.imshow(X, extent=extent, interpolation='bicubic',
                   vmin=0, vmax=1, **kwargs)
    return im


def gradient_bar(ax, x, y, width=0.5, bottom=0, cmap = cm.Blues_r):
    for left, top in zip(x, y):
        right = left + width
        gradient_image(ax, extent=(left, right, bottom, top),
                       cmap=cmap, cmap_range=(0.2,1.0),direction=1)

#plot(randn(5))
from matplotlib.colors import LinearSegmentedColormap
cdict = {'red':   [[0.0,  0.5, 0.5],
                   [0.5,  1.0, 1.0],
                   [1.0,  1.0, 1.0]],
         'green': [[0.0,  0.0, 0.0],
                   [0.5,  0.0, 0.0],
                   [1.0,  1.0, 1.0]],
         'blue':  [[0.0,  0.5, 0.5],
                   [0.5,  1.0, 1.0],
                   [1.0,  1.0, 1.0]]}
newcmp = LinearSegmentedColormap('Magentas', segmentdata=cdict, N=256)

plt.register_cmap(cmap=newcmp)

def plot_syninputs(fig,ax):
    xs = linspace(0,120,1000)
    ton, toff = (0.5,17.0)
    yglu = exp(-xs/toff)-exp(-xs/ton); yglu = yglu/max(yglu)
    ax.plot(concatenate(([-10],xs)),concatenate(([0],yglu)),'-',c=redcolor,linewidth= 2,alpha=0.5,label='NMDA')
    ax.set(autoscale_on=False)


    ton, toff = (0.1,1.8)
    trise = ton*toff/(ton-toff)*log(ton/toff)
    yglu = exp(-xs/toff)-exp(-xs/ton); yglu = yglu/max(yglu)
    ax.plot(concatenate(([-10],xs)),concatenate(([0],yglu))-1.015,'-',c=redcolor,linewidth= 2,label='AMPA')

    ton, toff = (0.5,15.0)
    trise = ton*toff/(ton-toff)*log(ton/toff)
    ygaba = exp(-xs/toff)-exp(-xs/ton); ygaba= ygaba/max(ygaba)
    ax.plot(concatenate(([-10],xs)),concatenate(([0],ygaba))-1.015*2,'g-',linewidth= 2,label='GABA')
    #gradient_bar(ax,[20],[1.25],3*toff,1.05,cmap = cm.Greens_r)
    ax.set_yticks(ticks=array([0,1.]))
    #gradient_bar(ax,[0],[1.27+.2],3*toff,1.27,cmap = newcmp)


    ax.set_aspect('auto')
    ax.set_xlim(-2,50)
    ax.set_ylim(-2.1,1.15)
    ax.set_xlabel('t (\si{\milli\second})')
    ax.set_ylabel('$\\textrm{G}/\\textrm{G}_{max}$',rotation=0, 
               horizontalalignment='left',
               verticalalignment='top')
    ax.yaxis.set_label_coords(0,1.05)
    #     ax.legend()   
    ax.annotate('NMDA',xy = (40,0.16))
    ax.annotate('AMPA',xy = (40,-0.95))
    ax.annotate('GABA',xy = (40,-1.88))


# In[11]:


#Jupyters: Proper Spine Model v2 - Main experiments

with open(folderstore+"All_baseline_datasetv3_"+condition+".pickle","rb") as f:
    vavgT,mesT,vtracsT,vtracsDT,vtracsST,CtracsT,me2T,_,dataT = pickle.load(f)

with open(folderstore+"All_baseline_depinhv3_"+condition+".pickle","rb") as f:
    messh,iPSDsh,posish,mesDiSI,spdataI = pickle.load(f)
   


# In[12]:


modeldt = 0.05
bands = 1


# In[13]:


im = Image.open('../Neuron_persp9.png')
height = im.size[1]
im = np.array(im).astype(np.float) / 255
fig = figure(figsize=(15,11.5))
gs = mpl.gridspec.GridSpec(3, 4,  wspace=0.35, hspace=.35) # 2x3 grid

ax0 = fig.add_subplot(gs[:1, 0:2]) # first full col
axl0 = fig.add_subplot(gs[1:, 0:2]) # first full col
ax1 = fig.add_subplot(gs[0, 2]) # first row, second col
ax3 = fig.add_subplot(gs[0, 3]) # first row, third col
ax2 = fig.add_subplot(gs[1, 2]) # 2nd row, second col
ax4 = fig.add_subplot(gs[1, 3]) # 2nd row, 3rd col
axl1 = fig.add_subplot(gs[2, 2]) # 3nd row, second col
axl2 = fig.add_subplot(gs[2, 3]) # 3nd row, 3rd col

plot_syninputs(fig,axl0)

ax0.imshow(im)
ax0.set_axis_off()

# Plots in fig EPSP
s1 = dataT['nPSD']==1.0
s2 = dataT['nPSD']==2.0

plot_trace(vtracsT[:,s1],arange(vtracsT.shape[0])*modeldt-100+5-25,ax1,c='C0',band = bands)
plot_trace(vtracsDT[:,s1],arange(vtracsDT.shape[0])*modeldt-100+5-25,ax1,c='C1',band = bands)
plot_trace(vtracsST[:,s1],arange(vtracsST.shape[0])*modeldt-100+5-25,ax1,c='C2',band = bands)
plot_trace(vtracsT[:,s2],arange(vtracsT.shape[0])*modeldt-100+5-25,ax1,c='r',band = 0,linestyle=':')





# Plot in fig Delta V - V_dend
s0 = dataT['A2']>0
s1 = ~s0
ax2.plot(abs(mesT[s1,0]),mesT[s1,1]/(mesT[s1,2]),'.',label="Spine head",alpha=0.5)
ax2.plot(abs(mesT[s0,0]),mesT[s0,1]/(mesT[s0,2]),'.',c=redcolor,label="Spine head",alpha=0.5)

# ax2.plot(abs(mesT[s1,0]),mesT[s1,1]-mesT[s1,2],'.',label="Spine head",alpha=0.5)
# ax2.plot(abs(mesT[s0,0]),mesT[s0,1]-mesT[s0,2],'.',c=redcolor,label="Spine head",alpha=0.5)


# Plot in fig Calcium
ax4.plot(abs(dataT['A1'][s1]),mesT[s1,4],'.',label="Spine head",alpha=0.5)
ax4.plot(abs(dataT['A1'][s0]),mesT[s0,4],'.',c=redcolor,label="Spine head",alpha=0.5)

# Plot in fig Delta V max
EL0 = -70
ax3.plot(abs(dataT['A1'][s1]),mesT[s1,1],'C0.',label="Spine head - SiS",alpha=0.5)
ax3.plot(abs(dataT['A1'][s0]),mesT[s0,1],'.',c=redcolor,label="Spine head - DiS",alpha=0.5)
ax3.plot(abs(dataT['A1']),me2T[:,3]-EL0,'C1.',label="Dendritic shaft",alpha=0.5)
ax3.plot(abs(dataT['A1']),mesT[:,3],'C2.',label="Soma",alpha=0.5)
#ax3.plot(abs(mesT[s0,-3]/1e-3),mesT[s0,3],'r.',label="Spine head",alpha=0.5)

#mVav = mVsshI[:,2].mean()
seldisi = spdataI['A2']>0
axl2.plot(posish,messh[:,3],'C2.',alpha=0.5,label='Shaft')#,label='Axo-dendritic')

axl2.plot(spdataI['Dss'][seldisi],mesDiSI[seldisi,3],'.',c=redcolor,alpha=0.5,label='DiS')#,label='Axo-spinous')


axl2.set_ylabel('$\\V_\\textrm{\Large max, soma}$\n(\si{\milli\\volt})',rotation=0, 
               horizontalalignment='left',
               verticalalignment='top')
axl2.set_xlabel('Distance to soma (\si{\micro\meter})')


#xt = linspace(min(mVsshI[:,-1]),max(mVsshI[:,-1]))
#axl2.plot(xt,interp(xt,xtshI,ytshI),'C2')

#axl2.plot(xt,interp(xt,xtspI,ytspI),'r')

#axl2.plot(xt,yt,'r',label='DiS')
#axl2.legend(title='Excitation',loc=(0.5,0.65))
#axl2.set_ylabel('$V_{\max}$ in soma',rotation = 0)

#axl2.fill_between(xt,yt+yt2,yt-yt2,color = 'r',alpha=0.7,band = 0)


s1v = dataT['nPSD']==1.0
s2v = dataT['nPSD']==2.0
axl1.plot(dataT['Dss'][s1v],mesT[s1v,3],'C0.',alpha=0.5,label='SiS')
axl1.plot(dataT['Dss'][s2v],mesT[s2v,3],'.',c=redcolor,alpha=0.5,label='DiS')

axl1.set_ylabel('$\\V_\\textrm{\Large max, soma}$\n(\si{\milli\\volt})',rotation=0, 
               horizontalalignment='left',
               verticalalignment='top')
axl1.set_xlabel('Distance to soma (\si{\micro\meter})')

xt = linspace(min(dataT['Dss']),max(dataT['Dss']))
#axl1.plot(xt,interp(xt,xtDiS,ytDiS),'r')
#axl1.plot(xt,interp(xt,xtSiS,ytSiS),'C0')#
#axl1.plot(xt,xt*0+getint(ytSiS)[0],'C0')

#axl2.set_ylabel('$V_{\max}$ in soma',rotation = 0)

#axl2.fill_between(xt,yt+yt2,yt-yt2,color = 'r',alpha=0.7,band = 0)


# Accesories

# scalebar(ax1,40,-65,5,2,xlab = '5 ms', ylab = '2 mV  ', color = 'k')
#ax1.set_axis_off()
ax2.set_xlabel("$\\R_\\textrm{\Large Neck}\, (\si{\Mohm})$")
# ax2.set_ylabel("$\\frac{\\V_\\textrm{\Large max, shaft}}{\\V_\\textrm{\Large max, spine}}$",rotation=0, 
#                horizontalalignment='left',
#                verticalalignment='top')
ax2.set_ylabel("$\\frac{\\V_\\textrm{\Large max, spine}}{\\V_\\textrm{\Large max, shaft}}$",rotation=0, 
               horizontalalignment='left',
               verticalalignment='top')

ax4.set_xlabel("$\\A_\\textrm{\Large ePSD}$ (\si{\square\micro\meter})")
ax4.set_ylabel("$[\\textrm{Ca}^{2+}]_\\textrm{\Large max}\, $(\si{\micro\Molar})",rotation=0, 
               horizontalalignment='left',
               verticalalignment='top')
#ax3.set_xlabel("- $I_{s,max}$(pA)")
ax3.set_xlabel("$\\A_\\textrm{\Large ePSD}$ (\si{\square\micro\meter})")

ax3.set_ylabel("$\Delta\\V_\\textrm{\Large max}\,$ (\si{\milli\\volt})",rotation=0, 
               horizontalalignment='left',
               verticalalignment='top')
ax3.tick_params(direction="in")
ax4.tick_params(which='both',direction="in")
ax2.tick_params(direction="in")
axl1.tick_params(direction="in")
axl2.tick_params(direction="in")
ax1.tick_params(which='both',direction="in")
axl0.tick_params(direction="in")

ax3.set_xticks(ticks=arange(0,1.5,0.5))
ax4.set_xticks(ticks=arange(0,1.5,0.5))

ax2.xaxis.set_minor_locator(FixedLocator(arange(100,600,100)))
ax2.xaxis.set_major_locator(FixedLocator(arange(0,600,200)))
ax2.xaxis.set_tick_params(which='minor',direction='in')
#ax4.xaxis.set_label_coords(0.5,-0.12)

#ax3.xaxis.set_label_coords(0.5,-0.12)
#ax2.xaxis.set_label_coords(0.5,-0.12)
#ax4.xaxis.set_label_coords(0.5,-0.12)

ax1.yaxis.set_label_coords(0.02,.97)
# ax2.yaxis.set_label_coords(0.02,.96)
ax2.yaxis.set_label_coords(0.06,.96)
ax3.yaxis.set_label_coords(0.02,.97)
ax4.yaxis.set_label_coords(0.035,.97)
axl2.yaxis.set_label_coords(0.02,.96)
axl1.yaxis.set_label_coords(0.02,.96)


#ax4.yaxis.set_label_coords(-.15,0.5)

ax3.legend(loc = (-2.55,0.45))
axl1.legend(title='EPSP',loc=(0.58,0.52))
axl2.legend(title='IPSP',loc=(0.52,0.52))

ax1.set_xlim(23-25,45-25)
ax1.set_ylim(-70.5,-54)
ax1.set_xlabel('t (\si{\milli\second})')
ax1.set_ylabel('V (\si{\milli\\volt})',rotation = 0 , 
               horizontalalignment='left',
               verticalalignment='top')


ax1.xaxis.set_minor_locator(FixedLocator([5,15]))
ax1.xaxis.set_major_locator(FixedLocator([0,10,20]))
ax1.yaxis.set_minor_locator(FixedLocator([-65,-55]))
ax1.yaxis.set_major_locator(FixedLocator([-70,-60,-50]))
ax1.xaxis.set_label_coords(0.5,-0.155)

extent = list(ax2.axis())
ax2.set_ylim(extent[2],extent[3]*1.17)
axl1.set_ylim(0,8.5)
axl2.set_ylim(0,25)

axl0.spines['top'].set_color('none')
axl0.spines['right'].set_color('none')

bboxor = ax0.figbox._points
dy =  bboxor[1,1]-bboxor[0,1]
dx = bboxor[1,0]-bboxor[0,0]
bboxn = [[bboxor[0,0]+dx*.45,bboxor[0,1]]-dy*.93,[bboxor[1,0]+dx*.15,bboxor[1,1]-dy*.21]]
ax0.set_position(Bbox(bboxn))


bboxor = axl0.figbox._points
dy =  bboxor[1,1]-bboxor[0,1]
dx = bboxor[1,0]-bboxor[0,0]
bboxn = [[bboxor[0,0]+dx*.15,bboxor[0,1]],[bboxor[1,0]-dx*.1,bboxor[1,1]-dy*.3]]
axl0.set_position(Bbox(bboxn))

ax2.set_ylim(1,4)
ax2.set_xlim(0,500)

#ax4.set_ylim(0,4.5)
axs = [ax0,ax1,ax3,ax2,ax4,axl1,axl2,axl0]
pos = zeros((len(axs),2))
pos[:,1] = 1-pos[:,1]
pos[:3,1] = 1.02
pos[0,1] = 1.127

pos[3:,1] = 1.15
pos[1:,0] = -0.08
pos[-1,:] = (-.052,1.1)
numbering_panels(axs,pos,labels=['A1','B','C','D','E','F','G','A2'])
fig.savefig(folderoutput+"Figure_Model_1v31.png",dpi = 300, tight_layout = True)

#avefig("Figure_Model_1bx2.pdf",dpi = 300, tight_layout = True)


# In[14]:


def gradient_image(ax, extent, direction=0.3, cmap_range=(0, 1), **kwargs):
    """
    Draw a gradient image based on a colormap.

    Parameters
    ----------
    ax : Axes
        The axes to draw on.
    extent
        The extent of the image as (xmin, xmax, ymin, ymax).
        By default, this is in Axes coordinates but may be
        changed using the *transform* kwarg.
    direction : float
        The direction of the gradient. This is a number in
        range 0 (=vertical) to 1 (=horizontal).
    cmap_range : float, float
        The fraction (cmin, cmax) of the colormap that should be
        used for the gradient, where the complete colormap is (0, 1).
    **kwargs
        Other parameters are passed on to `.Axes.imshow()`.
        In particular useful is *cmap*.
    """
    phi = direction * np.pi / 2
    v = np.array([np.cos(phi), np.sin(phi)])
    X = np.array([[v @ [1, 0], v @ [1, 1]],
                  [v @ [0, 0], v @ [0, 1]]])
    a, b = cmap_range
    X = a + (b - a) / X.max() * X
    im = ax.imshow(X, extent=extent, interpolation='bicubic',
                   vmin=0, vmax=1, **kwargs)
    return im


def gradient_bar(ax, x, y, width=0.5, bottom=0, cmap = cm.Blues_r):
    for left, top in zip(x, y):
        right = left + width
        gradient_image(ax, extent=(left, right, bottom, top),
                       cmap=cmap, cmap_range=(0.35,1.0),direction=1)

#plot(randn(5))
#ax = gca()
#ax.set(autoscale_on=False)
#gradient_bar(ax,[0],[0],3.,1.0)
#ax.set_aspect('auto')  


# In[15]:


from matplotlib.colors import LinearSegmentedColormap


# In[16]:


from matplotlib.colors import LinearSegmentedColormap
cdict = {'red':   [[0.0,  0.5, 0.5],
                   [0.5,  1.0, 1.0],
                   [1.0,  1.0, 1.0]],
         'green': [[0.0,  0.0, 0.0],
                   [0.5,  0.0, 0.0],
                   [1.0,  1.0, 1.0]],
         'blue':  [[0.0,  0.5, 0.5],
                   [0.5,  1.0, 1.0],
                   [1.0,  1.0, 1.0]]}
newcmp = LinearSegmentedColormap('Magentas', segmentdata=cdict, N=256)

plt.register_cmap(cmap=newcmp)

def plot_syninputs(fig,ax):
    xs = linspace(0,120,1000)
    ton, toff = (0.5,17.0)
    yglu = exp(-xs/toff)-exp(-xs/ton); yglu = yglu/max(yglu)
    ax.plot(concatenate(([-10],xs)),concatenate(([0],yglu)),'-',c=redcolor,linewidth= 2,alpha=0.5,label='NMDA')
    ax.set(autoscale_on=False)


    ton, toff = (0.1,1.8)
    trise = ton*toff/(ton-toff)*log(ton/toff)
    yglu = exp(-xs/toff)-exp(-xs/ton); yglu = yglu/max(yglu)
    ax.plot(concatenate(([-10],xs)),concatenate(([0],yglu))-1,'-',c=redcolor,linewidth= 2,label='AMPA')

    ton, toff = (0.5,15.0)
    trise = ton*toff/(ton-toff)*log(ton/toff)
    ygaba = exp(-xs/toff)-exp(-xs/ton); ygaba= ygaba/max(ygaba)
    ax.plot(concatenate(([-10],xs)),concatenate(([0],ygaba))-1-1,'g-',linewidth= 2,label='GABA')
    gradient_bar(ax,[20],[1.25],3*toff,1.05,cmap = cm.Greens_r)
    ax.set_yticks(ticks=array([0,1.]))
    gradient_bar(ax,[0],[1.27+.2],3*toff,1.27,cmap = newcmp)


    ax.set_aspect('auto')
    ax.set_xlim(-2,100)
    ax.set_ylim(-2.1,1.5)
    ax.set_xlabel('$t$ (ms)')
    ax.set_ylabel('$G/G_{max}$')
    ax.yaxis.set_label_coords(-0.03,0.4)
    ax.legend()
    ax.annotate('$\Delta t$',xy = (7,1.1))
    ax.annotate('',xy = (6,1.15),xytext = (0,1.15),arrowprops=dict(arrowstyle="<-"))
    ax.annotate('',xy = (14,1.15),xytext = (20,1.15),arrowprops=dict(arrowstyle="<-"))


# In[17]:


figure(figsize=(8,8))

xs = linspace(0,120,1000)

ton, toff = (0.5,17.0)
yglu = exp(-xs/toff)-exp(-xs/ton); yglu = yglu/max(yglu)


ton, toff = (0.1,1.8)
trise = ton*toff/(ton-toff)*log(ton/toff)
print('AMPA',trise)
yglu = exp(-xs/toff)-exp(-xs/ton); yglu = yglu/max(yglu)

ton, toff = (0.5,15.0)
trise = ton*toff/(ton-toff)*log(ton/toff)
print(trise)
ygaba = exp(-xs/toff)-exp(-xs/ton); ygaba= ygaba/max(ygaba)




xs = linspace(0,120,1000)
ton, toff = (0.5,17.0)
yglu = exp(-xs/toff)-exp(-xs/ton); yglu = yglu/max(yglu)


ton, toff = (0.1,1.8)
trise = ton*toff/(ton-toff)*log(ton/toff)
yglu = exp(-xs/toff)-exp(-xs/ton); yglu = yglu/max(yglu)

ton, toff = (0.5,15.0)
trise = ton*toff/(ton-toff)*log(ton/toff)
ygaba = exp(-xs/toff)-exp(-xs/ton); ygaba= ygaba/max(ygaba)


# In[18]:



def plot_syninputsdelta(fig,ax):
    xs = linspace(0,120,1000)

    ton, toff = (0.5,17.0)
    yglu = exp(-xs/toff)-exp(-xs/ton); yglu = yglu/max(yglu)
    ax.plot(concatenate(([-10],xs)),concatenate(([0],yglu)),'-',c=redcolor,linewidth= 2,alpha=0.5,label='NMDA')
    ax.set(autoscale_on=False)

    dy = 0.05
    gradient_bar(ax,[0],[1.25+.2-dy],2*toff,1.25-dy,cmap = newcmp)



    ton, toff = (0.1,1.8)
    trise = ton*toff/(ton-toff)*log(ton/toff)
    yglu = exp(-xs/toff)-exp(-xs/ton); yglu = yglu/max(yglu)
    ax.plot(concatenate(([-10],xs)),concatenate(([0],yglu))-1.015,'-',c=redcolor,linewidth= 2,label='AMPA')

    ton, toff = (0.5,15.0)
    trise = ton*toff/(ton-toff)*log(ton/toff)
    ygaba = exp(-xs/toff)-exp(-xs/ton); ygaba= ygaba/max(ygaba)
    ax.plot(concatenate(([-10],xs+20)),concatenate(([0],ygaba))-1.015*2,'g-',linewidth= 2,label='GABA')
    gradient_bar(ax,[20],[1.23-dy],2*toff,1.03-dy,cmap = cm.Greens_r)
    ax.set_yticks(ticks=array([0,1.]))



    ax.set_aspect('auto')

    ax.set_xlim(-2,70)
    ax.set_ylim(-2.1,1.6)
    ax.set_xlabel('t (\si{\milli\second})')
    ax.set_ylabel('$\\textrm{G}/\\textrm{G}_{max}$',rotation=0, 
               horizontalalignment='left',
               verticalalignment='top')
    ax.yaxis.set_label_coords(0,1.05)
    ax.vlines(20,-1.015*2,1.03-dy,linestyle='--')
    #ax.legend(loc = (0.6,.885))
    ax.annotate('$\Delta \\t$',xy = (6.5,1.05-dy))
    ax.annotate('NMDA',xy = (50,0.1))
    ax.annotate('AMPA',xy = (50,-0.95))
    ax.annotate('GABA',xy = (50,-1.85))
    ax.annotate('',xy = (6,1.13-dy),xytext = (0,1.13-dy),arrowprops=dict(arrowstyle="<-"))
    ax.annotate('',xy = (14,1.13-dy),xytext = (20,1.13-dy),arrowprops=dict(arrowstyle="<-"))


# In[21]:


#Jupyter notebook: Proper Spine Model - Inhibition.ipynb
data = {}

with open(folderstore+"dis_baselinev3_"+condition+".pickle","rb") as f:
    vavg,mes,vtracs,vtracsD,vtracsS,Ctracs = pickle.load(f)
data["dis_baseline"] = [vtracs,Ctracs]
with open(folderstore+"dis_ga_glu05v3_"+condition+".pickle","rb") as f:
    vavg,mes,vtracs,vtracsD,vtracsS,Ctracs = pickle.load(f)
data["dis_gaglu"] = [vtracs,Ctracs]
with open(folderstore+"dis_glu_ga05v3_"+condition+".pickle","rb") as f:
    vavg,mes,vtracs,vtracsD,vtracsS,Ctracs = pickle.load(f)
data["dis_gluga"] = [vtracs,Ctracs]
#with open("sis_ga_glu.pickle","rb") as f:
#    vavg,mes,vtracs,vtracsD,vtracsS,Ctracs = pickle.load(f)
data["sis_gaglu"] = [vtracs,Ctracs]
#with open("sis_glu_ga.pickle","rb") as f:
#    vavg,mes,vtracs,vtracsD,vtracsS,Ctracs = pickle.load(f)
data["sis_gluga"] = [vtracs,Ctracs]

with open(folderstore+"inhibition_v3_"+condition+".pickle","rb") as f:
    tdels,inhtimDis,_,tauDis,tauDis2,inhtimDism = pickle.load(f)

# Jupyter notebook: Proper Spine Model - Inhibition outside
with open(folderstore+"inhibition_v3_outPSD_"+condition+".pickle","rb") as f:
    tdels,inhtimSis,_,tauSis,tauSis2,inhtimSism = pickle.load(f)


# In[28]:


#height = im.size[1]
#im = np.array(im).astype(np.float) / 255
im = Image.open('../spine-inh_sketch_bothm.png')
height = im.size[1]
im = np.array(im).astype(np.float) / 255

fig = figure(figsize=(15,7.))
gs = mpl.gridspec.GridSpec(6,8,wspace=3, hspace=3) # 2x3 grid
ax0 = fig.add_subplot(gs[:2, :2]) # first full col
axl0 = fig.add_subplot(gs[2:, :2]) # first full col

ax1 = fig.add_subplot(gs[:3, 2:4]) # first row, second col
ax3 = fig.add_subplot(gs[3:, 2:4]) # 2nd row, second col

ax5 = fig.add_subplot(gs[:3, 4:6]) # first row, third col
ax0b = fig.add_subplot(gs[3:, 4:6]) # first full col
ax4 = fig.add_subplot(gs[3:, 6:]) # 2nd row, 3rd col
ax2 = fig.add_subplot(gs[:3, 6:]) # 2nd row, 3rd col

# Sketch
plot_syninputsdelta(fig,axl0)

ax0.imshow(im)
ax0.set_axis_off()

col_Inh = 'indigo'
col_ctr = redcolor

# Fig EPSP GABA-Glu dt = -5ms
vtracs = data["dis_baseline"][0]+70
plot_trace(vtracs,arange(vtracs.shape[0])*modeldt-100-20,ax1,c=col_ctr,band = bands)

vtracs = data["dis_gaglu"][0]+70
cG = col_Inh
plot_trace(vtracs,arange(vtracs.shape[0])*modeldt-100-20,ax1,c=cG,band = bands)
ax1.set_xlim(15-20,50-20)
#ax1.set_axis_off()

# Fig EPSP Glu-GABA dt = +5ms

vtracs = data["dis_baseline"][0]+70
#plot_trace(vtracs,arange(vtracs.shape[0])*.1-200+10,ax3,c='k')
#vtracs = data["dis_gaglu"][1]
plot_trace(vtracs,arange(vtracs.shape[0])*modeldt-100-20,ax3,c=col_ctr,band = bands,label='EPSP')

#vtracs = data["dis_baseline"][1]
#plot_trace(vtracs,arange(vtracs.shape[0])*.1-200+10,ax3,c='k')
vtracs = data["dis_gluga"][0]+70
plot_trace(vtracs,arange(vtracs.shape[0])*modeldt-100-20,ax3,c=cG,band = bands,label='EPSP \n\& IPSP')
ax3.set_xlim(15-20,50-20)
#ax1.set_ylim(-67,-20)
#ax3.set_axis_off()
#plot_trace(vtracsD,arange(vtracsD.shape[0])*.1-200+10,ax2,c='C1')
#plot_trace(vtracsS,arange(vtracsS.shape[0])*.1-200+10,ax2,c='C2')
#ax3.text(50,1,"$Ca^{2+}$")

vm,sl,sv = 1-inhtimSis[:,0],1-inhtimSis[:,1],1-inhtimSis[:,2]
ax2.plot(tdels,vm,'C0.-',label='Dendritic\n inhibition')
ax2.fill_between(tdels,sl,sv,alpha=0.5)

vm,sl,sv = 1-inhtimSis[:,3*5],1-inhtimSis[:,3*5+1],1-inhtimSis[:,3*5+2]
ax4.plot(tdels,vm,'C0.-',label='Dendritic\n inhibition')
ax4.fill_between(tdels,sl,sv,alpha=0.5,color='C0')
#ax2.set_axis_off()

vm,sl,sv = 1-inhtimDis[:,0],1-inhtimDis[:,1],1-inhtimDis[:,2]
ax2.plot(tdels,vm,'C1.-',label='Spinous\n inhibition')
ax2.fill_between(tdels,sl,sv,color='C1',alpha=0.5)

vm,sl,sv = 1-inhtimDis[:,3*5],1-inhtimDis[:,3*5+1],1-inhtimDis[:,3*5+2]
ax4.plot(tdels,vm,'C1.-',label='Spinous\n inhibition')
ax4.fill_between(tdels,sl,sv,color='C1',alpha=0.5)

vm,sl,sv = 1-inhtimSis[:,1*5],1-inhtimSis[:,1*5+1],1-inhtimSis[:,1*5+2]
ax5.plot(tdels,vm,'C0.-',label='Dendritic\n inhibition')
ax5.fill_between(tdels,sl,sv,alpha=0.5)

vm,sl,sv = 1-inhtimDis[:,1*5],1-inhtimDis[:,1*5+1],1-inhtimDis[:,1*5+2]
ax5.plot(tdels,vm,'C1.-',label='Spinous\n inhibition')
ax5.fill_between(tdels,sl,sv,color='C1',alpha=0.5)


#ax0btwin = ax0b.twinx()
#vm,sv = tauDis[:,2],tauDis[:,3]
#ax0btwin.plot(tdels,vm,'C1.-',label='In spine head')
##ax4twin.fill_between(tdels,vm+sv,vm-sv,color='C1',alpha=0.5)

#vm,sl,sv = tauSis[:,0],tauSis[:,1],tauSis[:,2]
#vm,sl,sv = tauSis2[:,0],tauSis2[:,1],tauSis2[:,2]
#vm,sl,sv = tauSis[:,5],tauSis[:,6],tauSis[:,7]
vm,sl,sv = tauSis2[:,15],tauSis2[:,16],tauSis2[:,17]
ax0b.plot(tdels,vm,'C0.-',label='Dendritic\n inhibition')
ax0b.fill_between(tdels,sl,sv,alpha=0.5)

#vm,sl,sv = tauDis[:,0],tauDis[:,1],tauDis[:,2]
#vm,sl,sv = tauDis2[:,0],tauDis2[:,1],tauDis2[:,2]
#vm,sl,sv = tauDis[:,5],tauDis[:,6],tauDis[:,7]
vm,sl,sv = tauDis2[:,15],tauDis2[:,16],tauDis2[:,17]
ax0b.plot(tdels,vm,'C1.-',label='Spinous\n inhibition')
ax0b.fill_between(tdels,sl,sv,color='C1',alpha=0.5)






ax2.spines['right'].set_color('none')
#ax2.spines['bottom'].set_color('none')
ax2.set_yticks(ticks=[0.0,0.1,0.2,0.3])
#ax2.set_xticks(ticks=arange(-50,70,20))
ax2.xaxis.set_minor_locator(FixedLocator(arange(-50,70,20)))
ax2.xaxis.set_major_locator(FixedLocator(arange(-40,60,20)))
ax2.tick_params(axis='x', pad = 10)
#ax2.yaxis.tick_right()
# Eliminate upper and right axes
ax2.spines['right'].set_position(('data',0.0))
ax2.spines['top'].set_position(('data',0.0))
#ax2.annotate('$V/V_{\\varnothing}$\n in spine',xy = (20,0.9))
ax2.spines['top'].set_color('none')
#ax2.xaxis.set_visible(False)
#ax2.set_ylim(.68,1.)


ax5.spines['right'].set_color('none')
#ax5.spines['bottom'].set_color('none')
ax5.spines['top'].set_color('none')
ax5.set_yticks(ticks=[0.0,0.1,0.20,0.3])
#ax2.set_xticks(ticks=arange(-50,70,20))
ax5.xaxis.set_minor_locator(FixedLocator(arange(-50,70,20)))
ax5.xaxis.set_major_locator(FixedLocator(arange(-40,60,20)))
ax5.tick_params(axis='x', pad = 10)
#ax5.yaxis.tick_right()
# Eliminate upper and right axes
#ax5.spines['right'].set_position(('data',0.0))
ax5.spines['top'].set_position(('data',0.0))
ax5.spines['top'].set_color('none')
#ax5.xaxis.set_visible(False)

ax4.spines['right'].set_color('none')
#ax4.spines['bottom'].set_color('none')

#ax4.xaxis.set_major_locator(FixedLocator(arange(-40,60,20)))

#ax4.yaxis.tick_right()
#ax4.xaxis.tick_top()

# Eliminate upper and right axes
#ax4.spines['right'].set_position(('data',0.0))
ax4.spines['bottom'].set_position(('data',0.00))
ax4.spines['top'].set_color('none')
#ax4.xaxis.set_visible(False)
#ax4.xaxis.tick_top()
ax4.set_ylabel('$S_{[{\\textrm Ca}^{2+}]}$',rotation = 0 , 
               horizontalalignment='left',
               verticalalignment='top')


ax0b.spines['top'].set_color('none')
ax0b.spines['right'].set_color('none')
#ax0b.spines['bottom'].set_color('none')
ax0b.yaxis.set_minor_locator(FixedLocator(arange(3,11,2)))
ax0b.yaxis.set_major_locator(FixedLocator(arange(4,12,2)))

#ax2.set_xticks(ticks=arange(-50,70,20))
#ax0b.tick_params(axis='x', pad = 10)
#ax0b.yaxis.tick_right()
#ax0b.xaxis.tick_top()
# Eliminate upper and right axes
#ax0b.spines['right'].set_position(('data',0.65))
#ax0b.spines['top'].set_color('none')
#ax0b.xaxis.set_visible(False)
ax0b.set_ylabel('$\\tau/\\tau_\infty$',rotation = 0 , 
               horizontalalignment='left',
               verticalalignment='top')


ax1.set(autoscale_on=False)
ax3.set(autoscale_on=False)
ax1.set(ylim=(-8,21),autoscale_on=False)
ax3.set(ylim=(-8,21),autoscale_on=False)
#ax1.plot(linspace(50,100),linspace(10,20)*0+2)

ddyg = 3
gradient_bar(ax1,[20-20],[-0.3-ddyg],17.*2,-2.5-ddyg,cmap = newcmp)
gradient_bar(ax1,[15-20],[-2.8-ddyg],15.*2,-5-ddyg,cmap = cm.Greens_r)

gradient_bar(ax3,[20-20],[-0.3-ddyg],17.*2,-2.5-ddyg,cmap = newcmp)
gradient_bar(ax3,[25-20],[-2.8-ddyg],15.*2,-5.0-ddyg,cmap = cm.Greens_r)

ax1.set_aspect('auto')
ax3.set_aspect('auto')
# scalebar(ax3,45,15,5,5,'$5$ ms','$5$ mV')
# scalebar(ax1,45,15,5,5,'$5$ ms','$5$ mV')

#gs.tight_layout(fig, rect=[.5,.5, 1, 1], h_pad=0.)#pad=0.4, w_pad=0.5, h_pad=4.0)
#print(xy1,sel)
#subplot_tool()

ax1.set_ylabel('$\\V$ (\si{\milli\\volt})',rotation = 0 , 
               horizontalalignment='left',
               verticalalignment='top')
ax3.set_ylabel('$\\V$ (\si{\milli\\volt})',rotation = 0 , 
               horizontalalignment='left',
               verticalalignment='top')
ax1.set_xlabel('$\\t$ (\si{\milli\second}) ')
ax3.set_xlabel('$\\t$ (\si{\milli\second}) ')



ax0b.set_xlabel('$\Delta \\t$ (\si{\milli\second}) ')
ax4.set_xlabel('$\Delta \\t$ (\si{\milli\second}) ')
ax5.set_xlabel('$\Delta \\t$ (\si{\milli\second}) ')
ax2.set_ylabel('$S_\\V$ in spine',rotation = 0 , 
               horizontalalignment='left',
               verticalalignment='top')
ax5.set_ylabel('$S_\\V$ in dend',rotation = 0 , 
               horizontalalignment='left',
               verticalalignment='top')

ax4.yaxis.set_label_coords(0.015,1.0)
ax0b.yaxis.set_label_coords(0.035,1.03)
ax2.yaxis.set_label_coords(0.035,1.)
ax5.yaxis.set_label_coords(0.035,1.)
ax1.yaxis.set_label_coords(0.035,1.2)
ax3.yaxis.set_label_coords(0.035,1.2)

ax4.xaxis.set_minor_locator(FixedLocator(arange(-50,70,20)))
ax4.xaxis.set_major_locator(FixedLocator(arange(-40,60,20)))
ax0b.xaxis.set_minor_locator(FixedLocator(arange(-50,70,20)))
ax0b.xaxis.set_major_locator(FixedLocator(arange(-40,60,20)))
ax0b.yaxis.set_major_locator(FixedLocator(arange(.6,1.2,.2)))
##ax0b.yaxis.set_minor_locator(FixedLocator(arange(4,14,4)))
ax0b.yaxis.set_minor_locator(FixedLocator(arange(.5,1.1,0.2)))

#
#ax2.set_xticklabels([])
#ax5.set_xticklabels([])
ax5.xaxis.set_label_coords(0.5,-0.15)
ax2.set_xlabel('$\Delta t$ (\si{\milli\second}) ')
ax2.xaxis.set_label_coords(0.5,-0.15)
ax2.tick_params(which='minor',direction='in')
ax2.tick_params(which='major',direction='in')
ax4.tick_params(which='minor',direction='in')
ax4.tick_params(which='major',direction='in')
ax5.tick_params(which='minor',direction='in')
ax5.tick_params(which='major',direction='in')
ax0b.tick_params(which='minor',direction='in')
ax0b.tick_params(which='major',direction='in')

ax3.tick_params(which='minor',direction='in')
ax3.tick_params(which='major',direction='in')
ax1.tick_params(which='minor',direction='in')
ax1.tick_params(which='major',direction='in')
ax2.tick_params(pad=2.5)
ax5.tick_params(pad=2.5)

ax0b.set_ylim(.0,1.05)

xl0,xlf = -50,20
ax0b.set_xlim(xl0,xlf)
ax5.set_xlim(xl0,xlf)
ax4.set_xlim(xl0,xlf)
ax2.set_xlim(xl0,xlf)

ax2.set_ylim(0.,.8)
ax5.set_ylim(0.,.8)
# ax1.annotate('$\Delta t=-5$ ms',xy = (30,7))
# ax3.annotate('$\Delta t=+5$ ms',xy = (35,7))
#ax0b.legend(loc = (-0.47,0.05))
ax2.legend(loc = (-0.95,0.28))

ax3.legend(loc=(0.45,1.88))

ax3.annotate("", xy=(2+1,20-2), xytext=(2+11,20-2),
    arrowprops=dict(arrowstyle="<->"))

#ax3.arrow(50+22,.9*exp(-1.0),60,0.,head_width=0.2,head_length=20)

#ax3.vlines(20+3,0.1,18,linestyle='--')
ax3.vlines(2+1,9.76,19,linestyle='--')
ax3.vlines(2+11,2.44,19,linestyle='--')

ax3.annotate('$\\tau$',xy = (27-20,21-2))

ax1.annotate('DiS: $V_{\\rm spine}$',xy = (40,35))

ax0b.set_ylim(0.0,1.05)


xst = linspace(0,200,2000)
biexp = exp(-xs/10.0)-exp(-xs/1.3)
biexp = biexp/max(biexp)
#ax1.plot(xs+20,biexp*25,'--',linewidth=2)


# bboxor = ax3.figbox._points
# dy = bboxor[0,1]*0.48
# bboxn = [[bboxor[0,0],bboxor[0,1]-dy],[bboxor[1,0],bboxor[1,1]-dy*.5]]
# ax3.set_position(Bbox(bboxn))

# bboxor = ax1.figbox._points
# bboxn = [[bboxor[0,0],bboxor[0,1]-dy],[bboxor[1,0],bboxor[1,1]-dy*.5]]
# ax1.set_position(Bbox(bboxn))
bboxor = ax3.figbox._points
dy = bboxor[1,1]-bboxor[0,1]
fac = 0.275
bboxn = [[bboxor[0,0],bboxor[0,1]-dy*fac],[bboxor[1,0],bboxor[1,1]-dy*fac]]
ax3.set_position(Bbox(bboxn))

bboxor = ax1.figbox._points
fac = 0.28
bboxn = [[bboxor[0,0],bboxor[0,1]-dy*fac],[bboxor[1,0],bboxor[1,1]-dy*fac]]
ax1.set_position(Bbox(bboxn))

bboxor = ax0.figbox._points
dy =  bboxor[1,1]-bboxor[0,1]
dx = bboxor[1,0]-bboxor[0,0]
bboxn = [[bboxor[0,0]+dx*0.15,bboxor[0,1]]-dy*.2,[bboxor[1,0]+dx*0.15,bboxor[1,1]+dy*.15]]
# bboxn = [[bboxor[0,0],bboxor[0,1]-dy],[bboxor[1,0],bboxor[1,1]-dy]]
ax0.set_position(Bbox(bboxn))


axl0.spines['top'].set_color('none')
axl0.spines['right'].set_color('none')
bboxor = axl0.figbox._points
dy =  bboxor[1,1]-bboxor[0,1]
dx = bboxor[1,0]-bboxor[0,0]
bboxn = [[bboxor[0,0],bboxor[0,1]],[bboxor[1,0]+dx*.25,bboxor[1,1]]]
axl0.set_position(Bbox(bboxn))


bboxor = ax5.figbox._points
dy =  bboxor[1,1]-bboxor[0,1]
dx = bboxor[1,0]-bboxor[0,0]
bboxn = [[bboxor[0,0]-dx*.15,bboxor[0,1]],[bboxor[1,0]-dx*.15,bboxor[1,1]]]
ax5.set_position(Bbox(bboxn))

bboxor = ax0b.figbox._points
dy =  bboxor[1,1]-bboxor[0,1]
dx = bboxor[1,0]-bboxor[0,0]
bboxn = [[bboxor[0,0]-dx*.15,bboxor[0,1]],[bboxor[1,0]-dx*.15,bboxor[1,1]]]
ax0b.set_position(Bbox(bboxn))

# ax3.tick_params(axis='x', which='major', pad=27)
# ax1.tick_params(axis='x', which='major', pad=27)


ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax1.spines['bottom'].set_position(('data',0.00))
ax3.spines['bottom'].set_position(('data',0.00))
ax3.spines['left'].set_bounds(0.0,24)
ax1.spines['left'].set_bounds(0.0,24)

ax5.vlines(0,0.0,0.35,linestyle='--')
ax2.vlines(0,0.0,0.35,linestyle='--')
ax0b.vlines(0,0.0,1.2,linestyle='--')
ax4.vlines(0,0.0,0.2,linestyle='--')

ax1.xaxis.set_label_coords(0.95,0.455)
ax3.xaxis.set_label_coords(0.95,0.455)

ax4.set_yticks(ticks=[0.0,0.1, 0.2])
ax4.set_ylim(0.,0.3)
#ax4.yaxis.set_minor_locator(FixedLocator(arange(0.05,0.55,0.1)))

axs = [ax0,ax1,ax5,ax2,ax0b,ax4,axl0,ax3]
pos = zeros((len(axs),2))
pos[:,1] = 1-pos[:,1]
pos[1:,1] = 1.10
pos[1:,0] = -0.08
pos[4:,1] = 1.15

pos[0,0] = 0.11
pos[0,1] = 1.22
pos[1,1] = 1.38

pos[-2,0] = 0.0
pos[-1,0] = 0.02
pos[1,0] = 0.02
pos[-1,1] = 1.3

numbering_panels(axs,pos,labels=['A1','B1','C','D','E','F','A2','B2'])


fig.savefig(folderoutput+'f2v3.png',dpi = 300,tight_layout = True)

#savefig('f2x2.pdf',dpi = 300,tight_layout = True)
im = Image.open('../Neuron_persp9.png')
height = im.size[1]
im = np.array(im).astype(np.float) / 255
fig = figure(figsize=(15,11.5))
gs = mpl.gridspec.GridSpec(3, 4,  wspace=0.35, hspace=.35) # 2x3 grid

ax0 = fig.add_subplot(gs[:1, 0:2]) # first full col
axl0 = fig.add_subplot(gs[1:, 0:2]) # first full col
ax1 = fig.add_subplot(gs[0, 2]) # first row, second col
ax3 = fig.add_subplot(gs[0, 3]) # first row, third col
ax2 = fig.add_subplot(gs[1, 2]) # 2nd row, second col
ax4 = fig.add_subplot(gs[1, 3]) # 2nd row, 3rd col
axl1 = fig.add_subplot(gs[2, 2]) # 3nd row, second col
axl2 = fig.add_subplot(gs[2, 3]) # 3nd row, 3rd col

plot_syninputs(fig,axl0)

ax0.imshow(im)
ax0.set_axis_off()

# Plots in fig EPSP
s1 = dataT['nPSD']==1.0
s2 = dataT['nPSD']==2.0





# Plot in fig Delta V - V_dend
s0 = dataT['A2']>0
s1 = ~s0
ax2.plot(abs(mesT[s1,0]),(mesT[s1,2])/mesT[s1,1],'.',label="Spine head",alpha=0.5)
ax2.plot(abs(mesT[s0,0]),(mesT[s0,2])/mesT[s0,1],'.',c=redcolor,label="Spine head",alpha=0.5)

# ax2.plot(abs(mesT[s1,0]),mesT[s1,1]-mesT[s1,2],'.',label="Spine head",alpha=0.5)
# ax2.plot(abs(mesT[s0,0]),mesT[s0,1]-mesT[s0,2],'.',c=redcolor,label="Spine head",alpha=0.5)


# Plot in fig Calcium
ax4.plot(abs(dataT['A1'][s1]),mesT[s1,4],'.',label="Spine head",alpha=0.5)
ax4.plot(abs(dataT['A1'][s0]),mesT[s0,4],'.',c=redcolor,label="Spine head",alpha=0.5)

# Plot in fig Delta V max
EL0 = -70
ax3.plot(abs(dataT['A1'][s1]),mesT[s1,1],'C0.',label="Spine head - SiS",alpha=0.5)
ax3.plot(abs(dataT['A1'][s0]),mesT[s0,1],'.',c=redcolor,label="Spine head - DiS",alpha=0.5)
ax3.plot(abs(dataT['A1']),me2T[:,3]-EL0,'C1.',label="Dendritic shaft",alpha=0.5)
ax3.plot(abs(dataT['A1']),mesT[:,3],'C2.',label="Soma",alpha=0.5)
#ax3.plot(abs(mesT[s0,-3]/1e-3),mesT[s0,3],'r.',label="Spine head",alpha=0.5)

#mVav = mVsshI[:,2].mean()
seldisi = spdataI['A2']>0
s1v = dataT['nPSD']==1.0
s2v = dataT['nPSD']==2.0

xt = linspace(min(dataT['Dss']),max(dataT['Dss']))

# Accesories


ax2.set_xlabel("$\\R_\\textrm{\Large Neck}\, (\si{\Mohm})$")
# ax2.set_ylabel("$\\frac{\\V_\\textrm{\Large max, shaft}}{\\V_\\textrm{\Large max, spine}}$",rotation=0, 
#                horizontalalignment='left',
#                verticalalignment='top')
ax2.set_ylabel("$\\frac{\\V_\\textrm{\Large max, shaft}}{\\V_\\textrm{\Large max, spine}}$",rotation=0, 
               horizontalalignment='left',
               verticalalignment='top')

ax4.set_xlabel("$\\A_\\textrm{\Large ePSD}$ (\si{\square\micro\meter})")
ax4.set_ylabel("$[\\textrm{Ca}^{2+}]_\\textrm{\Large max}\, $(\si{\micro\Molar})",rotation=0, 
               horizontalalignment='left',
               verticalalignment='top')
#ax3.set_xlabel("- $I_{s,max}$(pA)")
ax3.set_xlabel("$\\A_\\textrm{\Large ePSD}$ (\si{\square\micro\meter})")

ax3.set_ylabel("$\Delta\\V_\\textrm{\Large max}\,$ (\si{\milli\\volt})",rotation=0, 
               horizontalalignment='left',
               verticalalignment='top')
ax3.tick_params(direction="in")
ax4.tick_params(which='both',direction="in")
ax2.tick_params(direction="in")
axl0.tick_params(direction="in")

ax3.set_xticks(ticks=arange(0,1.5,0.5))
ax4.set_xticks(ticks=arange(0,1.5,0.5))

ax2.xaxis.set_minor_locator(FixedLocator(arange(100,600,100)))
ax2.xaxis.set_major_locator(FixedLocator(arange(0,600,200)))
ax2.xaxis.set_tick_params(which='minor',direction='in')
#ax4.xaxis.set_label_coords(0.5,-0.12)

#ax3.xaxis.set_label_coords(0.5,-0.12)
#ax2.xaxis.set_label_coords(0.5,-0.12)
#ax4.xaxis.set_label_coords(0.5,-0.12)

# ax2.yaxis.set_label_coords(0.02,.96)
ax2.yaxis.set_label_coords(0.06,.96)
ax3.yaxis.set_label_coords(0.02,.97)
ax4.yaxis.set_label_coords(0.035,.97)
axl2.yaxis.set_label_coords(0.02,.96)
axl1.yaxis.set_label_coords(0.02,.96)


#ax4.yaxis.set_label_coords(-.15,0.5)

ax3.legend(loc = (-2.55,0.45))


extent = list(ax2.axis())
ax2.set_ylim(extent[2],extent[3]*1.17)

axl0.spines['top'].set_color('none')
axl0.spines['right'].set_color('none')

bboxor = ax0.figbox._points
dy =  bboxor[1,1]-bboxor[0,1]
dx = bboxor[1,0]-bboxor[0,0]
bboxn = [[bboxor[0,0]+dx*.45,bboxor[0,1]]-dy*.93,[bboxor[1,0]+dx*.15,bboxor[1,1]-dy*.21]]
ax0.set_position(Bbox(bboxn))


bboxor = axl0.figbox._points
dy =  bboxor[1,1]-bboxor[0,1]
dx = bboxor[1,0]-bboxor[0,0]
bboxn = [[bboxor[0,0]+dx*.15,bboxor[0,1]],[bboxor[1,0]-dx*.1,bboxor[1,1]-dy*.3]]
axl0.set_position(Bbox(bboxn))

ax2.set_ylim(0,1.05)
ax2.set_xlim(0,500)



axt = axl1

plot_trace(vtracsT[:,s1],arange(vtracsT.shape[0])*modeldt-100+5-25,axt,c='C0',band = bands)
#plot_trace(vtracsDT[:,s1],arange(vtracsDT.shape[0])*modeldt-100+5-25,axt,c='C1',band = bands)
#plot_trace(vtracsST[:,s1],arange(vtracsST.shape[0])*modeldt-100+5-25,axt,c='C2',band = bands)
axt.tick_params(which='both',direction="in")
axt.yaxis.set_label_coords(0.02,.97)
axt.set_xlim(23-25,45-25)
axt.set_ylim(-70.5,-60)
axt.set_xlabel('t (\si{\milli\second})')
axt.set_ylabel('V (\si{\milli\\volt})',rotation = 0 , 
               horizontalalignment='left',
               verticalalignment='top')
axt.xaxis.set_minor_locator(FixedLocator([5,15]))
axt.xaxis.set_major_locator(FixedLocator([0,10,20]))
axt.yaxis.set_minor_locator(FixedLocator([-65,-55]))
axt.yaxis.set_major_locator(FixedLocator([-70,-60,-50]))
axt.xaxis.set_label_coords(0.5,-0.155)


axt = axl2

#plot_trace(vtracsT[:,s1],arange(vtracsT.shape[0])*modeldt-100+5-25,axt,c='C0',band = bands)
plot_trace(vtracsDT[:,s1],arange(vtracsDT.shape[0])*modeldt-100+5-25,axt,c='C1',band = bands)
#plot_trace(vtracsST[:,s1],arange(vtracsST.shape[0])*modeldt-100+5-25,axt,c='C2',band = bands)
axt.tick_params(which='both',direction="in")
axt.yaxis.set_label_coords(0.02,.97)
axt.set_xlim(23-25,45-25)
axt.set_ylim(-70.5,-60)
axt.set_xlabel('t (\si{\milli\second})')
axt.set_ylabel('V (\si{\milli\\volt})',rotation = 0 , 
               horizontalalignment='left',
               verticalalignment='top')
axt.xaxis.set_minor_locator(FixedLocator([5,15]))
axt.xaxis.set_major_locator(FixedLocator([0,10,20]))
axt.yaxis.set_minor_locator(FixedLocator([-65,-55]))
axt.yaxis.set_major_locator(FixedLocator([-70,-60,-50]))
axt.xaxis.set_label_coords(0.5,-0.155)


axt = ax1

#plot_trace(vtracsT[:,s1],arange(vtracsT.shape[0])*modeldt-100+5-25,axt,c='C0',band = bands)
#plot_trace(vtracsDT[:,s1],arange(vtracsDT.shape[0])*modeldt-100+5-25,axt,c='C1',band = bands)
plot_trace(vtracsST[:,s1],arange(vtracsST.shape[0])*modeldt-100+5-25,axt,c='C2',band = bands)
axt.tick_params(which='both',direction="in")
axt.yaxis.set_label_coords(0.02,.97)
axt.set_xlim(23-25,45-25)
axt.set_ylim(-70.5,-60)
axt.set_xlabel('t (\si{\milli\second})')
axt.set_ylabel('V (\si{\milli\\volt})',rotation = 0 , 
               horizontalalignment='left',
               verticalalignment='top')
axt.xaxis.set_minor_locator(FixedLocator([5,15]))
axt.xaxis.set_major_locator(FixedLocator([0,10,20]))
axt.yaxis.set_minor_locator(FixedLocator([-65,-55]))
axt.yaxis.set_major_locator(FixedLocator([-70,-60,-50]))
axt.xaxis.set_label_coords(0.5,-0.155)

ax2.yaxis.set_label_coords(0.05,1.05)

#ax4.set_ylim(0,4.5)
axs = [ax0,ax1,ax3,ax2,ax4,axl1,axl2,axl0]
for axt in axs:
    try:
        axt.spines['top'].set_visible(False)
        axt.spines['right'].set_visible(False)
    except:
        pass
    
pos = zeros((len(axs),2))
pos[:,1] = 1-pos[:,1]
pos[:3,1] = 1.02
pos[0,1] = 1.127

pos[3:,1] = 1.15
pos[1:,0] = -0.08
pos[-1,:] = (-.052,1.1)
numbering_panels(axs,pos,labels=['A1','B','C','D','E','F','G','A2'])
savefig("Figure_Model_1v31splitting.png",dpi = 300, tight_layout = True)
fig.savefig(folderoutput+'f1v31splitting.png',dpi = 300,tight_layout = True)
