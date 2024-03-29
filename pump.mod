TITLE Sodium-potassium pump
 

UNITS {
        (molar) = (1/liter)
        (pA) =  (picoamp)
	(mV) =	(millivolt)
        (uS) =  (micromho)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {v FROM -100 TO 50 WITH 50 (mV)}

NEURON {
	SUFFIX pump
	USEION na READ nai WRITE ina
	USEION k WRITE ik
	RANGE  inapump,ipumpmax,n,km,ina,ik
 
}


PARAMETER {
        dt  (ms)
        nai (mM)
        celsius = 35 (degC)
        ipumpmax  = 0.0092 (mA/cm2)
        km = 10.0 (mM)
        n  = 1.5
        
        
}

ASSIGNED { 
           ina	   (mA/cm2)
           ik	   (mA/cm2)
           inapump (mA/cm2)
}

BREAKPOINT {
    inapump = ipumpmax*(1/(1 + pow(km/nai,n)))
	:ina = 3.0*inapump
	:ik = -2.0*inapump
	ina = 3.0*0.00021
	ik = -2.0*0.00021
}
