TITLE ERG channel from Tucker et al 2012

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(S)  =  (siemens)
}

NEURON {
	SUFFIX erg
	USEION k WRITE ik
	RANGE gergbar,ik,vshift,gk
}

PARAMETER {
	ek = -90.0 (mV)
	gergbar= 10e-6	(S/cm2) :  must be explicitly def. in hoc
	vshift = 0 (mV)
}

STATE { 
	ierg <1e-10> oerg <1e-10>  cerg <1e-10>
}

ASSIGNED {
	v	(mV)
	ik	(mA/cm2)
	alpha_a (1/ms)
	beta_a (1/ms)
	alpha_i (1/ms)
	beta_i (1/ms)
	gk (S/cm2)
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	gk  = gergbar*oerg
	ik = gk*(v-ek)
}



KINETIC kin {
	rates(v-vshift)
	~ oerg <-> ierg (alpha_i,beta_i)
	~ oerg <-> cerg (beta_a,alpha_a)
	
	CONSERVE oerg + cerg + ierg = 1

}

INITIAL {
	rates(v-vshift)
	SOLVE kin
	STEADYSTATE sparse
	:cerg=0.9
	:oerg=0.01
	:ierg=0.09
}



FUNCTION exponential(v(mV),a(1/ms),b(1/mV)) {
	exponential = a*exp(b*v)
}

UNITSOFF

PROCEDURE rates(v) {
	:alpha_a = exponential(v,0.0036,0.0759)
	:beta_a = exponential(v,1.2523e-05,-0.0671)
	:alpha_i = exponential(v,91.11,0.1189)
	:beta_i = exponential(v,12.6,0.0733)
	
	alpha_a = exponential(v,0.0061,0.1055)
	beta_a = exponential(v,2.1469e-5,-0.057)
	alpha_i = exponential(v,227.775,0.1123)
	beta_i = exponential(v,21.0,0.0712)
}

UNITSON
