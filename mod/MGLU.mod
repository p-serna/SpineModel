TITLE    MGLU synapse for nucleus accumbens model
: see comments below

NEURON {
	POINT_PROCESS MGLU
	:SUFFIX MGLU_er
	RANGE tau, beta, gamma, scale, ip3itemp, t1, ip3i, ical, ip3min, spkcnt
}

UNITS {
	(molar) = (1/liter)	
	(mM)	= (millimolar)
	(nA) = (nanoamp)
    (mA) = (milliampere)

	
    (um) = (micron)
    FARADAY = (faraday) (coulomb)
    PI = (pi) (1)
}

PARAMETER {
	gamma = 5e-6 (mM/ms2)
	tau = 220 (ms)
	beta = 0.2 (1/ms)
	ip3min = 0.0001 (mM)
	ip3i0=0.00005 (mM)
}


ASSIGNED {
	t1 (ms)
	spkcnt
	scale			: scale allows the current to be scaled by weight

	ip3i (mM)
}



STATE{
	ip3alpha (mM/ms)
	ip3itemp (mM)
}


INITIAL {
	t1 = -10000
	spkcnt=0
	ip3itemp=ip3i0
	ip3alpha=0
	scale = 1
}

BREAKPOINT {
  	SOLVE states METHOD derivimplicit
	:cnexp  euler  derivimplicit
}

DERIVATIVE states {
	ip3itemp' =  scale*gamma*(t-t1) * exp(-(t-t1)/tau)  - beta * (ip3itemp - ip3min)
	:ip3itemp' =  scale*ip3alpha  - beta * (ip3itemp - ip3min)
	ip3i=ip3itemp
}

NET_RECEIVE(weight) {
	:ip3alpha = ip3alpha +ip3alpha1
	:ip3alpha=ip3alpha*gamma*(t-t1) * exp(-(t-t1)/tau)
	
	: store the spike time
	t1 = t

	spkcnt = spkcnt +1
	scale = weight
}
