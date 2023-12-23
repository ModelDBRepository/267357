
from neuron import h, gui
import numpy as np

# because of the way neuron objects are not destroyed easily, 
# only set one variable to 1 for any given run
RUN=0
NULLCLINE=0
BIFURCATION=0
NANULL=0
BALANCED=1
BALARRAY=0

prefix = 'balanced_state_1006'
h.load_file("template_simple.hoc") # single compartment template
h.load_file("fixnseg.hoc")
verbose= 0 
v_init = -62.6

p = dict({})	# keep parameters in a dictionary allos for passing
#""" #Figure 5 parameters			# arbitrary numbers/order of parameters to class creation
p['na_cond'] =  550.0e-6
p['kdr_cond'] = 460.0e-6
p['ca_cond'] = 5.196e-6
p['kca_cond'] = 30.0e-6 # figure 5
p['a_cond_s'] = 570.0e-6
p['a_cond_p'] = 300e-6#285.0e-6
p['a_cond_d'] = 266.0e-6
p['erg_cond'] = 100e-6
p['ca_pump'] = 0.0012 
p['kleak'] = 0.5e-6 #compensated by erg 'leak' # used for bifurcation/nullclines as fixed value k+ conductance

p['iamp']=0


p['nmdaamp'] = 10.0e-7
p['gabaamp'] = 20.0e-7
#"""
""" # Figure 4
p['na_cond'] =  500.0e-6
p['kdr_cond'] = 445.0e-6 # at 440 for 4c, 445 for 4b
p['ca_cond'] = 5.196e-6
p['kca_cond'] = 60.0e-6 
p['a_cond_s'] = 570.0e-6
p['a_cond_p'] = 285.0e-6#285.0e-6
p['a_cond_d'] = 266.0e-6
p['erg_cond'] = 60e-6 #50
p['ca_pump'] = 0.0012 #0.00191
p['kleak'] = 0.5e-6 #compensated by erg 'leak' # used for bifurcation/nullclines as fixed value k+ conductance

p['iamp']=0

# not used by 4
p['nmdaamp'] = 12.5e-7
p['gabaamp'] = 12.5e-7
#"""

g_celsius = 35
h.celsius = g_celsius
import new_geomnseg as geom # python version of geomnseg - call after 
							# conductances assigned for multi compartmental models


class dopa:
	def __init__(self,i, params):
		SINGLE=1
		#h.pop_section()
		self.nrn = h.SIMPLE()

		self.nrn.soma.nseg = 1
		self.nrn.soma.diam = 2
		self.nrn.soma.L=10

	
		for s in self.nrn.all:
			#s.insert('nabalan')
			s.insert('hh3')
			s.gnabar_hh3 = params['na_cond']
			s.gkhhbar_hh3 = params['kdr_cond']
			#if s in self.nrn.proxDend:
			s.gkabar_hh3 = params['a_cond_p']
			#elif s in self.nrn.distDend:
				#s.gkabar_hh3 = a_cond_d
			#else:
				#s.gkabar_hh3 = a_cond_s	

			s.insert('pump')
			s.insert('leak')
			s.gkbar_leak = params['kleak']
			s.insert('cabalan')
			s.icapumpmax_cabalan = params['ca_pump']
			s.insert('cachan')
			s.mhalf_cachan = -45.0
			s.mslope_cachan = 5.0
			s.gcalbar_cachan = params['ca_cond']
			#s.gkbar_cachan = 0.0
			s.insert('kca')
			s.gkbar_kca = params['kca_cond']
			s.insert('erg')
			s.gergbar_erg = params['erg_cond']
			
		self.iapp = h.IClamp(0.5,sec=self.nrn.soma)
		self.iapp.dur = 1e9
		self.iapp.amp = params['iamp']
		self.iapp.delay = 0
		self.nc = h.NetCon(self.nrn.soma(0.5)._ref_v,None,sec=self.nrn.soma) 
		# creating a netcon forces a sufficiently small dt to resolve spikes when using parallelcontext
		self.nc.threshold=-20
		self.srec = h.Vector()
		self.nc.record(self.srec)

		#self.nc2 = h.NetCon(self.nrn.soma(0.5)._ref_v,None,sec=self.nrn.soma)
		#self.nc2.threshold=-60
		#self.srec2 = h.Vector()
		#self.nc2.record(self.srec2)
		"""
		self.i2 = h.IClamp(0.5,sec=self.nrn.soma)
		self.i2.dur = 100
		self.i2.amp = -params['iamp']
		self.i2.delay = 15000
		#"""
	def __del__(self):
			pass

def runfunc():

	cells = []
	
	p0=p.copy()
	
	cells.append(dopa(0,p0))
	
	p1 = p.copy()
	p1['kca_cond'] *= 0
	#p1['erg_cond'] *=1.25
	cells.append(dopa(0,p1))

	p2 = p.copy()
	p2['kca_cond'] *= 0
	#p2['erg_cond'] *=1.25
	p2['na_cond'] = 0
	
	cells.append(dopa(0,p2))

	
	# reproduction of oscillations from Canavier et al 2007
	vvecs = []
	ivecs=[]
	svecs = []
	i=0
	for c in cells:
		vvecs.append(h.Vector())
		ivecs.append(h.Vector())
		svecs.append(h.Vector())
		svecs.append(h.Vector())
		svecs.append(h.Vector())
		vvecs[i].record(c.nrn.soma(0.5)._ref_v,0.1)
		ivecs[i].record(c.nrn.soma(0.5)._ref_gk_erg,0.1)
		svecs[3*i].record(c.nrn.soma(0.5)._ref_oerg_erg,0.1)
		svecs[3*i+1].record(c.nrn.soma(0.5)._ref_ierg_erg,0.1)
		svecs[3*i+2].record(c.nrn.soma(0.5)._ref_h_hh3,0.1)
		i+=1

	cv=h.CVode()
	cv.active(1)
	#pc.nthread(6)
	#pc.set_maxstep(0.1)

	#h.finitialize()
	h.finitialize(v_init)
	h.t = 0
	h.dt = 0.1
	#print cells[0].nrn.soma.ena
	while h.t < 30000:
		cv.solve(h.t+1000)
		#h.fadvance()
		#if h.t%1000 < h.dt:
		print( h.t )
	
	print 'done single cell'
	fp = open('%sergtraces.dat'% prefix,'w')
	for i in range(len(vvecs[0])):
		if i < 20000:
			continue
		fp.write('%e  ' % (i*0.1))
		for j in range(len(cells)):
			fp.write('%e  %e  %e  %e  %e  ' %(vvecs[j].x[i], ivecs[j].x[i],svecs[3*j+2].x[i],svecs[3*j].x[i]*p1['erg_cond'],svecs[3*j].x[i]+svecs[3*j+1].x[i]))
		fp.write('\n')
			
	fp.close()
	return 0

if RUN:
	runfunc() # putting cell creation inside a function destroys them on completion

# create nullclines
cells = None

def nullclinefunc():
	ptemp = p.copy() # dictionaries cannot be copied using a = b
		 # eg ptemp = p is equivalent to to ptemp = *p in c
	fixedergcond = ptemp['erg_cond']
	ptemp['erg_cond'] = 0
	ptemp['kca_cond'] = 0
	#ptemp['na_cond'] = 0
	#ptemp['a_cond_p'] = 0

	baseline = ptemp['kleak'] # normal conductance
	condsweep = np.arange(baseline,50.1e-6,0.1e-6) # only uses values that produce points within vb
	ptemp['kleak'] = condsweep[0]
	# replacing slow potassium conductances with constant values
	# to produce bifurcation, nullclines.
	import nullcline_funcs as nf # this creates new cells
	vb=[-75,-20]
	# create erg nullcline
	v = -100
	#"""
	cell = dopa(0,ptemp)
	
	fp = open('%serg_steady.dat' % (prefix),'w')
	vdata = []
	cdata = []
	odata = []
	idata = []

	while v < 40.1:
		h.finitialize(v)
		c = cell.nrn.soma
		fp.write('%e  %e  %e  %e   %e\n' %(v, c.oerg_erg, c.cerg_erg, c.ierg_erg,c.oerg_erg+c.ierg_erg))
		vdata.append(v)
		cdata.append(c.cerg_erg)
		odata.append(c.oerg_erg)
		idata.append(c.ierg_erg)
		v+=1
	fp.close()

	from scipy.interpolate import interp1d
	cfunc = interp1d(vdata,cdata,kind='cubic')
	ofunc = interp1d(vdata,odata,kind='cubic')
	ifunc = interp1d(vdata,idata,kind='cubic')
	
	def oiratio(v):
		vtemp = v
		if v <= -100:
			vtemp=-100
		if v >= 40:
			vtemp=40
		
		return ofunc(vtemp)/ifunc(vtemp)
	
	area = h.area(0.5,sec=cell.nrn.soma)
	
	del cell
	#"""
	#return 0 # uncomment for

	gvpairs = []
	cells = []
	has_started = 0
	#balh = dopa(0,ptemp)
	for conductance in condsweep:
		if cells ==[]:
			test, cells = nf.v_root(dopa, cells=[],args=(0,ptemp),volt_bounds=vb,voltstep=1)
			#print 'fine'
		else:

			ptemp['kleak'] = conductance
			for c in cells:
				for s in c.nrn.all:
					s.gkbar_leak = conductance

			# quasi-newton root finder allows for finding root(s) in spite of negative slope regions
			test, cells = nf.v_root(dopa, cells=cells, args=(0,ptemp),volt_bounds=vb,voltstep=1) 
				
		if len(test) > 0 and not has_started:
			has_started = 1
		if len(test) < 1 and has_started:
			break
		print(conductance, test)
		
		for things in test:
			gvpairs.append([things,conductance-baseline])
			
	print('built v of gk')
	gvpairs.sort()
	fname = '%sslowk_v_nullcline' %prefix
	if ptemp['a_cond_p'] == 0:
		fname += '_noKv4'
	if ptemp['na_cond']==0:
		fname += '_noNa'

	fname += '.dat'
	fp = open(fname,'w')
	

	pc = h.ParallelContext()
	pc.nthread(1) 	# forces single thread computations required for the brute force jacobian
					# don't worry its actually fast unless your neuron has HUGE numbers of states
	
	
	
	
	ptemp = p.copy()
	ptemp['kca_cond'] = 0
	
	
	for pairs in gvpairs:
		max_eig, comp, ncomp, eigs = nf.stability_check(dopa,pairs[0],args=(0,ptemp),offset=-ptemp['iamp'], rundur=5000)
		print(pairs[0], max_eig, comp)
		oerg = pairs[1]/fixedergcond
		slow_erg = oerg+oerg/oiratio(pairs[0])
		
		fp.write('%e  %e  %e  %e  %e  %e\n' % (pairs[1],pairs[0],oerg, slow_erg, max_eig, comp)) # voltage, current, stability, complex component
	
	fp.close()

	del pc
	
	return 0

if NULLCLINE:
	nullclinefunc()
	
	



def bifurcationfunc():
	STEPLEN = 20000 # time in ms to reach equilibrium for IV curve, 
	tol=1 # how close extrema can be considered non-unique (in mV)
	
	import nullcline_funcs as nf

	iarray = np.arange(-3e-4, 3e-3, 1.0e-5)
	pc = h.ParallelContext()
	pc.nthread(64)
	cells = []
	garray = np.arange(0,10e-6,1e-7)
	
	ptemp = p.copy()
	ptemp['erg_cond'] = 0
	#ptemp['kca_cond'] = 0
	ptemp['kleak'] = 0.5e-6+garray[0]
	
	vb = [-80,-20]
	#f, varray, curarray, cells = nf.IV_curve(dopa,cells=[],args=(0,ptemp),volt_bounds=vb,rundur=STEPLEN,voltstep=0.1)
	
	ncells = len(cells) # we might as well keep the generated cells


	stim = []
	kicks = []
	vecs = []
	tvecs=[]
	i = 0
	h.t=0
	h.dt= 0.1
	#gp = open('demo.dat','w') # created for debugging 
	"""
	for current in iarray:
		ptemp['iamp'] = current
		if i >= ncells:
			cells.append(dopa(0,ptemp))
			ncells+=1
		else:
			cells[i].iapp.amp = current
		#"""
	for cond in garray:
		ptemp['kleak'] = cond+0.5e-6
		if i >= ncells:
			cells.append(dopa(0,ptemp))
			ncells+=1
		else:
			for s in cells[i].nrn.all:
				s.gkbar_leak = 0.5e-6+cond		
		#vecs.append(h.Vector())
		#tvecs.append(h.Vector())
		#vecs[i].record(cells[i].nrn.soma(0.5)._ref_v, 0.1)
		#tvecs[i].record(h._ref_t)
		
		# 'kick' knocks cell out of bistable state, reveals potenial oscillation
		kicks.append(h.IClamp(cells[i].nrn.soma(0.5)))
		kicks[i].amp = 10e-3
		kicks[i].dur = 100
		kicks[i].delay = STEPLEN/2.0
		
		#stim.append(h.IClamp(cells[i].nrn.soma(0.5)))
		#stim[i].amp = current
		#if current > 0:
		#	kicks[i].amp -= current # kick to 0 current
		#stim[i].delay = 0
		#stim[i].dur = 1e9
		i+=1
		

	vecs = []
	ivecs=[]
	svecs = []

	while len(cells) > len(garray):
		cells[-1] = None
		cells.pop()
		ncells-=1
		
		
	for c in cells:
		vecs.append(h.Vector())
		svecs.append(h.Vector())
		vecs[-1].record(c.nrn.soma(0.5)._ref_v,0.1)
		svecs[-1].record(c.nrn.soma(0.5)._ref_oerg_erg,0.1)

			

	#cv = h.CVode()
	#cv.active(1)
	h.finitialize(-70)
	print('integrating models')
	
	
	
	while h.t < STEPLEN:
		#cv.solve(h.t+100)
		pc.psolve(h.t+1000)
		print(h.t,cells[0].nrn.soma.v)
	
	#free up some un-needed currents
	
	del kicks

	
	fp = open('%serg_envelopes.dat' % prefix, 'w')
	
	print(len(cells), ncells, len(vecs)) # all should be same
	
	from scipy.signal import argrelextrema
	
	for i in range(ncells):
		#print len(vecs[i]),
		vecs[i].remove(0,10*(STEPLEN-5000)) # remove all but last 5 seconds
		#print len(vecs[i]), cells[i].iapp.i, vecs[i].x[0],vecs[i].x[len(vecs[i])-1],vecs[i].max(),vecs[i].min()
		pyvec = np.zeros(len(vecs[i]))
		pyvec = vecs[i].to_python(pyvec)
				
		"""if abs(iarray[i]) < 1e-7:
			for j in range(len(vecs[i])):
				gp.write('%e\n' % vecs[i].x[j])
			gp.close()
			quit()"""


		# note that due to how parallel context works, finding local extrema is unreliable
			
		print (garray[i],pyvec.max(),pyvec.min())

		fp.write('%e  %e\n' % (garray[i],pyvec.min()))
		if abs(pyvec.min()-pyvec.max()) > tol:
			fp.write('%e  %e\n' % (garray[i],pyvec.max()))
		
	fp.close()
	quit()
	del stim
	del vecs
	del tvecs
	ptemp['iamp'] = 0

	#quit()
	pc.nthread(1) # required for next bit as calls to cvode.f require single thread
	#"""
	
	fp = open('%sIV_erg.dat' %prefix , 'w')
	try:	
		del pc
		del cells
	except:
		pass
	cv = h.CVode()
	j=0
	for v in varray:
		try:
			i = f(v)
		except:
			i=curarray[j]
		
		ptemp['iapp'] = i
		max_eig, comp, ncomp, eigs = nf.stability_check(dopa,v ,args=(0,ptemp),offset=0,rundur=30000)
		print( float(i), float(v), max_eig, comp, ncomp )
		#print len(eigs)
		fp.write('%e  %e  %e  %e  %d\n' % (i,v,max_eig,comp,ncomp))
		j+=1
	fp.close()
	
	return 0
	
	
def gen_netstim(func,target,syn_type='gabaa',nspikes=1,duration=1e4,funcargs={},pstable=p,input_spikes=None,ncells=1,rtype='nexp'):
	rtype = 'nexp'
	if rtype not in ['poisson','nexp']:
		rtype = 'poisson'
		if VERBOSE:
			print 'type not supported, reverting to poisson'
	safety_interval = 1000
	#stim = h.NetStim()
	#print pstable['gabaamp']
	

	
	if target == None: # creates the spike train, but doesn't create a netcon, synapse
		synapse = None
	else:
		if syn_type not in ['nmda','gabab']:
			synapse = h.Exp2Syn(target(0.5))
			synapse.tau1=1
			synapse.tau2=6
			if syn_type=='gabaa':
				synapse.e = -70
			if syn_type=='ampa':
				synapse.e = 0
				synapse.tau2=3
		if syn_type == 'nmda':
			synapse = h.Exp2NMDA(target(0.5))
		if syn_type == 'gabab':
			#synapse = h.GABAb(target(0.5)) # this doesn't have right syntax
			synapse = h.Exp2Syn(target(0.5))
			synapse.tau1=50
			synapse.tau2=1000
			synapse.e = -90
	
	#print synapse
	if target == None:	
		ncobj = None
	else:	
		ncobj=h.NetCon(None, synapse)
		ncobj.delay=2
		if syn_type=='gabaa':
			ncobj.weight[0] = pstable['gabaamp']
		if syn_type=='nmda':
			ncobj.weight[0] = pstable['nmdaamp']
		if syn_type=='gabab':
			ncobj.weight[0] = pstable['gababamp']
		#print syn_type, ncobj.weight[0]		

	
	if input_spikes != None:
		if isinstance(input_spikes,list): # if it is a python array, just copy it
			spikes = input_spikes
			if target != None:
				return ncobj, spikes, synapse
			else:
				return spikes
		else: # if it is array like (eg numpy array or neuron vector)
			try:
				spikes = []

				for things in input_spikes:
					spikes.append(things)
				if target != None:
					return ncobj, spikes, synapse
				else:
					return spikes
			except:
				print 'input spikes expected to be an array object, generating from scratch'
				pass
		
	
	
	keys = funcargs.keys()
	execstr = 'freq = func(t'
	for key in keys:
	   execstr+= ',%s=%f' % (key, funcargs[key])
	execstr+=')'
	spikes = []
	sftinv = 1.0/safety_interval
	old=0
	for j in range(ncells):
		t=0
		while t < duration: # this can cause problems if the minimum frequency is very low.
			nspikes=1
			exec(execstr) # freq= func(t,args)
			#print freq
			if freq < sftinv:
				t+=10
				continue
			if rtype=='poisson':
				isi = np.random.poisson(1000.0/freq,nspikes)
			if rtype=='nexp':
				isi = np.random.exponential(1000.0/freq,nspikes)
			#print isi
			if max(isi) > safety_interval: # 
				t+=safety_interval
				continue
			cs = np.cumsum(isi)
			for things in cs:
				if things < safety_interval:
					spikes.append(things+t)
					old = things
				else:
					break
			t = t + old
	if ncells > 1:
		spikes.sort()
	
	print len(spikes), syn_type, funcargs
	#test1 = h.Vector(spikes).deriv(1)
	#for i in range(len(test1)):
		#print spikes[i],test1[i]
	if target != None:
		return ncobj, spikes, synapse
	else:
		return spikes
	
	
def stepfunc(t, basal=0,delay=0,duration=1000,pulseamp=0):
	if delay <= t <= delay+duration:
		return pulseamp
	else:
		return basal


def balanced():
	rec_list = []
	
	dt =0.1
	
	#p['kca_cond']=0*40e-6
	#p['erg_cond']=60e-6
	cell = dopa(0,p)

	
	TSTOP=40000
	
	cell.nc.record(cell.srec) # this might be doubling up on record commands

	rec_vectors = [h.Vector(),h.Vector()]
	i=2
	rec_vectors[0].record(cell.iapp,h._ref_t,dt,sec=cell.nrn.soma)
	rec_vectors[1].record(cell.iapp,cell.nrn.soma(0.5)._ref_v,dt,sec=cell.nrn.soma)
	#print rec_list
	for things in rec_list:
		try:
			rec_vectors.append(h.Vector())
			exec('rec_vectors[%d].record(blah,cell.nrn.soma(0.5)._ref_%s,dt,sec=cell.nrn.soma)' % (i,things) ) in locals()
			i+=1
		except:
			print 'failed to record %s'  % things
	

	cv = h.CVode()
	cv.active(1)
	
	spikes_gaba = gen_netstim(stepfunc, None,nspikes=10,duration=TSTOP,funcargs={'basal':1000,'delay':3000,'duration':1000,'pulseamp':0},ncells=1,syn_type='gabaa')
	spikes_nmda = gen_netstim(stepfunc, None,nspikes=10,duration=TSTOP,funcargs={'basal':100,'delay':7000,'duration':1000,'pulseamp':0},ncells=1,syn_type='nmda')
	
	
	gaba_syn = h.Exp2Syn(cell.nrn.soma(0.5))

	gaba_syn.tau1=1
	gaba_syn.tau2=6
	gaba_syn.e = -70
	
	nmda_syn = h.Exp2NMDA(cell.nrn.soma(0.5))
	
	
	gaba_nc = h.NetCon(None, gaba_syn)
	nmda_nc = h.NetCon(None, nmda_syn)
	
	gaba_nc.weight[0] = p['gabaamp']
	nmda_nc.weight[0] = p['nmdaamp']
	
	h.finitialize()
	
	# events need to come after finitialize
	for spikes in spikes_gaba:
		gaba_nc.event(spikes)
	for spikes in spikes_nmda:
		nmda_nc.event(spikes)
		
		
	cv.solve(TSTOP)
	
	sz = len(rec_vectors[0])
	fp = open('%sbalanced_state_%e_%e.dat' %(prefix,p['gabaamp'],p['nmdaamp']),'w')
	for i in range(sz):
		for vecs in rec_vectors:
			fp.write('%e  ' % vecs.x[i])
		fp.write('\n')
		
	fp.close()
	
	
		
	# create a single 


def balanced_array():
	rec_list = []
	
	dt =0.1
	
	#p['kca_cond']=40e-6
	#p['erg_cond'] = 40e-6
	#cell = dopa(0,p)

	nmda_array = np.arange(0,20.01e-7,0.5e-7)
	gaba_array = np.arange(0,40.01e-7,1e-7)
	ncells = len(nmda_array)*len(gaba_array)
	# 900 cells
	TSTOP=20000
	print ncells

	

	pc=h.ParallelContext()
	pc.nthread(64)
	
	spikes_gaba = gen_netstim(stepfunc, None,nspikes=10,duration=TSTOP,funcargs={'basal':1000,'delay':3000,'duration':0,'pulseamp':0},ncells=1,syn_type='gabaa')
	spikes_nmda = gen_netstim(stepfunc, None,nspikes=10,duration=TSTOP,funcargs={'basal':100,'delay':7000,'duration':0,'pulseamp':0},ncells=1,syn_type='nmda')
	# with no variation this can be done with netstim
	cells = []
	gaba_syns = []
	gaba_ncs = []
	nmda_syns= []
	nmda_ncs = []
	subobvs = []
	subnc = []
	subrec=[]
	i=0
	for nmda_cond in nmda_array:
		for gaba_cond in gaba_array:
			cells.append(dopa(0,p))
			subobvs.append(h.SecNet(cells[i].nrn.soma(0.5)))
			subobvs[i].thresh = -60
			subrec.append(h.Vector())
			subnc.append(h.NetCon(subobvs[i],None,sec=cells[i].nrn.soma))
			
			#subnc[i].threshold = 0 # checks for start of depolarizing wave
			subnc[i].record(subrec[i])
			
			
			
			cells[i].nc.record(cells[i].srec) # this might be doubling up on record commands
			gaba_syns.append(h.Exp2Syn(cells[i].nrn.soma(0.5)))
			gaba_syns[i].tau1=1
			gaba_syns[i].tau2=6
			gaba_syns[i].e = -70		
			nmda_syns.append(h.Exp2NMDA(cells[i].nrn.soma(0.5)))	
			gaba_ncs.append(h.NetCon(None, gaba_syns[i]))
			nmda_ncs.append(h.NetCon(None, nmda_syns[i]))
			gaba_ncs[i].weight[0] = gaba_cond
			nmda_ncs[i].weight[0] = nmda_cond
			#print(i)
			i+=1
	h.finitialize()
	
	# events need to come after finitialize
	for i in range(ncells):
		for spikes in spikes_gaba:
			gaba_ncs[i].event(spikes)
		for spikes in spikes_nmda:
			nmda_ncs[i].event(spikes)
	
	while h.t <TSTOP: # this can take a while without large numbers of cores	
		pc.psolve(h.t+100)
		print(h.t)
	
	isitemp = h.Vector()
	fp = open('%sbalarraydata.dat'%prefix,'w')
	for i in range(ncells):
		cells[i].srec.where('>',1000)# chop off initial transient
		subrec[i].where('>',1000)
		nspikes = len(cells[i].srec)
		nwaves = len(subrec[i])
		kill=0
		if nspikes>2:
			isitemp.deriv(cells[i].srec,1)
			freq = 1000.0/isitemp.mean()
			var = isitemp.var()/isitemp.mean()
			minisi = isitemp.min()
			isburst = 0
			burstiness = 0
			for spikes in isitemp:
				if spikes < 80:
					isburst=1
				if isburst and spikes > 150:
					isburst=0
				if isburst:
					burstiness+=1
			burstiness /= (1.0*nspikes)
			kill=0
		else:
			burstiness = 0
			freq = 0
			var = 0
			minisi=100
			kill=1
			
		print(nspikes,nwaves)
		if nwaves > 2 and nspikes > 2:
			temp=h.Vector()
			inburst = 0
			for j in range(nwaves-1):
				temp.where(cells[i].srec,'()',subrec[i].x[j],subrec[i].x[j+1])
				if len(temp) > 1:
					inburst += len(temp)
				if len(temp) == 0:
					nwaves -= 1
		else:
			inburst = nspikes
			nwaves = 1
		if nspikes > TSTOP/5000:
			inburst = inburst/float(nspikes)
		else:
			nspikes = 0
			inburst = 0

		#print(nspikes,nwaves)
		
		toggle = 0
		#nspikes SHOULD be > nwaves
		if 0.5<=(nspikes)/float(nwaves) <= 1.5:
			toggle = 1
		
		elif (nspikes)/float(nwaves) > 1.5:
			toggle = 2
		else:
			toggle = 0
		
		fp.write('%e  %e  %e  %e  %e  %e  %e  %e  %d\n' %(gaba_ncs[i].weight[0],nmda_ncs[i].weight[0],freq,var,burstiness,minisi,nspikes/float(nwaves),inburst,toggle))
		print('%e  %e  %e  %e  %e  %e  %e  %e  %d\n' %(gaba_ncs[i].weight[0],nmda_ncs[i].weight[0],freq,var,burstiness,minisi,nspikes/float(nwaves),inburst,toggle))

		#print('%e  %e  %e' %(gaba_ncs[i].weight[0],nmda_ncs[i].weight[0],len(cells[i].srec)/(0.001*TSTOP))   )
	fp.close()

	
if BIFURCATION:
	bifurcationfunc()
	

if BALANCED:
	balanced()
	
if BALARRAY:
	balanced_array()
quit()
