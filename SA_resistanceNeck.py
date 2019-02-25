%pylab
x = linspace(0,1,100)

xkcd()
figure(figsize = (8,6))

plot((x-0.4)/0.4*100,100*(1./(1-x**2)-1.0))
ylim(0,300)
xlabel("Increase in diameter of SA (%)")
ylabel("Increase Resistance of spine neck (%)")
savefig("SA_percentage_increaseRes.pdf")


figure(figsize = (8,6))
D0 = 1.0
DSA = 1.2
plot((x-0.4)/0.4*100,(D0/DSA)**2./(1-x**2))
ylim(.2,3)
xlabel("Increase in diameter of SA (%)")
ylabel("R$_{Neck}$ with SA /R$_{Neck}$ without SA  ")
xp = array([0.4,sqrt(1-(D0/DSA)**2),0.8])
plot((xp-0.4)/0.4*100,(D0/DSA)**2./(1-xp**2),'ko')
xp = linspace(-100,(xp[1]*.95-0.4)/0.4*100,10)
plot(xp,xp*0+1.0,'k--')
xlim(-100,150)

annotate(
    'If spine apparatus occupies 80% of\n the cross section of the spine neck,\n resistance is roughly 2 times\n the average resistance ',
    xy=(96, 1.96), arrowprops=dict(arrowstyle='-'), xytext=(-60, 2.3))

annotate(
    'It has to occupy ~55% to\n recover average resistance\n of spines without SA',
    xy=(38., 1.05), arrowprops=dict(arrowstyle='-'), xytext=(-60, 1.5))
    
annotate(
    'If SA blocks\n the spine neck,\n R$_{Neck}$ with SA$\\rightarrow\infty$',
    xy=(148., 2.98), arrowprops=dict(arrowstyle='-'), xytext=(50, 0.3))
    
savefig("SA_ResvsnonSA.pdf")


figure(figsize = (8,6))
D0 = 1.0
DSA = 1.2
plot(x*100,(D0/DSA)**2./(1-x**2))
ylim(.2,3)
xlabel("Diameter of SA (% of diameter spine neck)")
ylabel("R$_{Neck}$ with SA /R$_{Neck}$ without SA  ")
xp = array([0.4,sqrt(1-(D0/DSA)**2),0.8])
plot(xp*100,(D0/DSA)**2./(1-xp**2),'ko')
xp = linspace(-100,xp[1]*98,10)
plot(xp,xp*0+1.0,'k--')
xlim(0,100)

annotate(
    'If spine apparatus occupies 80% of\n the cross section of the spine neck,\n resistance is roughly 2 times\n the average resistance ',
    xy=((96/250+.4)*100, 1.96), arrowprops=dict(arrowstyle='-'), xytext=((-60/250+.4)*100, 2.3))

annotate(
    'It has to occupy ~55% to\n recover average resistance\n of spines without SA',
    xy=((38./250+.4)*100, 1.05), arrowprops=dict(arrowstyle='-'), xytext=((-60/250+.4)*100, 1.5))
    
annotate(
    'If SA blocks\n the spine neck,\n R$_{Neck}$ with SA$\\rightarrow\infty$',
    xy=((148./250+.4)*100, 2.98), arrowprops=dict(arrowstyle='-'), xytext=((50/250+.4)*100, 0.3))

#SBF-EM    
annotate(
    'Estimated diameter\n from SEM images',
    xy=(39.5, 0.785), arrowprops=dict(arrowstyle='-'), xytext=(20, 0.35))
    
savefig("SA_ResvsnonSAb.pdf")

