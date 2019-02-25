
COMMENT
T-type Ca channel 
ca.mod for ca current inside spines from Ly et al 2015
Uses fixed eca instead of GHK eqn
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX itL
	USEION ca READ eca WRITE ica
	RANGE m, h, gca, gbar
	RANGE minf, hinf, mtau, htau, inactF, actF
	GLOBAL  vshift,vmin,vmax, v12m, v12h, vwm, vwh, am, ah, vm1, vm2, vh1, vh2, wm1, wm2, wh1, wh2, ah2, ah3, ahh, ahh2, ahh3, vhh1, vhh2, whh1, whh2, vth
}

PARAMETER {
	gbar = 0.25974 (mho/cm2)	: 0.12 mho/cm2
	vshift = 0	(mV)		: voltage shift (affects all)

	cao  = 2.5	(mM)	        : external ca concentration
	cai		(mM)
						 
	v 		(mV)
	dt		(ms)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)

	v12m=47         	(mV)
	v12h=80         	(mV)
	vwm =5.0         	(mV)
	vwh=6.0         	(mV)
	am=1.84         	(mV)
	ah=7.6         	(mV)
    ah2=177         (mV)
    ah3=134        (mV)
	vm1=27         	(mV)
	vm2=72         	(mV)
	vh1=56.6         	(mV)
	vh2=99.5         	(mV)
	wm1=8         	(mV)
	wm2=20         	(mV)
	wh1=6.3         	(mV)
	wh2=5.6         	(mV)
    ahh = 3.10      (ms)
    ahh2= 3.683     (ms)
    ahh3= 46.34     (ms)
    vhh1= 37.9      (mV)
    vhh2= 70.6      (mV)
    whh1= 4.7       (mV)
    whh2= 8.6       (mV)
    vth = -60       (mV)
    


}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
	PI	= (pi) (1)
} 

ASSIGNED {
	ica 		(mA/cm2)
	gca		(pS/um2)
	eca		(mV)
	minf 		hinf
	mtau (ms)	htau (ms)
	tadj
}
 

STATE { m h }

INITIAL { 
	trates(v+vshift)
	m = minf
	h = hinf
}

BREAKPOINT {
        SOLVE states
        gca = gbar*m*m*m*h
	ica = gca * (v - eca)
} 

LOCAL mexp, hexp

PROCEDURE states() {
        trates(v+vshift)      
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
	VERBATIM
	return 0;
	ENDVERBATIM
}


PROCEDURE trates(v) {  
                      
        LOCAL tinc
        TABLE minf, mexp, hinf, hexp
	DEPEND dt	
	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable == 1

        tinc = -dt 

        mexp = 1 - exp(tinc/mtau)
        hexp = 1 - exp(tinc/htau)
}


PROCEDURE rates(v_) {  
        LOCAL  a, b

	minf = 1.0 / ( 1 + exp(-(v_+v12m)/vwm) )
	hinf = 1.0 / ( 1 + exp((v_+v12h)/vwh) )

    mtau = 1.0 / ( am  / (1.0+ exp(-(v_+vm1)/wm1) + exp(-(v_+vm2)/wm2) ) )
    
    htau = 1.0 / ( ah + ah2 / (1.0+ exp(-(v_+vh1)/wh1)) + ah3/(1.0+ exp(-(v_+vh2)/wh2) ) )
    if (v > vth) {
	htau = ( ahh + ahh2 / (1.0+ exp(-(v_+vhh1)/whh1))+ahh3 / (1.0 + exp((v_+vhh2)/whh2) ) ) 
	}
    

}



