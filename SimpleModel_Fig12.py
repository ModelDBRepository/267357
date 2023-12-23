
from neuron import h, gui
import numpy as np

RUN=1
NULLCLINE=1
BALANCED=0
BALARRAY=0

prefix = 'control_current_'

h.load_file("template_simple.hoc")
timecon = 0
back= 0 
sakmann=0
h.load_file("fixnseg.hoc")
verbose= 0 
v_init = -62.6
vcseries = 0
clamp = 0                     #switch for voltage clamp*/
restart = 0 #1*/                      # switch for initializing from state.old */
tstart = 0
if vcseries:
	tstop = 800 
else:
	tstop = 5000    #time in msec*/
if clamp and not vcseries:
	tstop= 600
if timecon:
	tstop = 800


nsyn = 45			# The number of synapses */
Dtmax = 1.0  
Dt = 1.00
if timecon:
	Dtmax = 1.0 
dt = 5e-1 #5e-4*/  
if timecon :
	dt = 0.01
nainit=  4.075
vsolder=v_init
vsold=v_init#PARAMETERS*/
na_cond =  550.0e-6
kdr_cond = 650.0e-6
ca_cond = 11.196e-6
kca_cond = 59.0e-6
a_cond_s = 570.0e-6
a_cond_p = 285.0e-6
a_cond_d = 266.0e-6

p = dict({})	# keep parameters in a dictionary allos for passing
				# arbitrary numbers/order of parameters to class creation
p['na_cond'] =  550.0e-6
p['kdr_cond'] = 550.0e-6 # at ~460 there is a transition into depblock
p['ca_cond'] =11.196e-6
p['kca_cond'] = 59.0e-6
p['a_cond_s'] = 570.0e-6
p['a_cond_p'] = 285e-6#285.0e-6
p['a_cond_d'] = 266.0e-6

p['vinit'] = -55



p['nmdaamp'] = 1.5e-6
p['gabaamp'] = 3e-5

p['iamp']=0
#stronger gA *1.28 =  729.6, 364.8, 340.48*/
iapl = 0 #in nA, -0.180nA=-180pA*/

global_ra = 40
global_cm = 1.0

g_celsius = 35
h.celsius = g_celsius
import new_geomnseg as geom # python version

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
			#s.gkbar_leak = params['kleak']
			s.insert('cabalan')
			s.fCa_cabalan=0.02
			s.icapumpmax_cabalan *=1.35
			s.insert('cachan')
			s.mhalf_cachan = -40.0
			s.mslope_cachan = 7.0
			s.gcalbar_cachan = params['ca_cond']
			#s.gkbar_cachan = 0.0
			s.insert('kca')
			s.gkbar_kca = params['kca_cond']
			#s.insert('erg')
			#s.gergbar_erg = params['erg_cond']
			
		self.iapp = h.IClamp(0.5,sec=self.nrn.soma)
		self.iapp.dur = 1e9
		self.iapp.amp = params['iamp']
		self.iapp.delay = 0
		self.nc = h.NetCon(self.nrn.soma(0.5)._ref_v,None,sec=self.nrn.soma) 
		# creating a netcon forces a sufficiently small dt to resolve spikes when using parallelcontext
		self.nc.threshold=-20
		self.srec = h.Vector()
		self.nc.record(self.srec)
		"""
		self.i2 = h.IClamp(0.5,sec=self.nrn.soma)
		self.i2.dur = 100
		self.i2.amp = -params['iamp']
		self.i2.delay = 15000
		#"""
	def __del__(self):
			pass
					
				
		#geom.geom_nseg(self.nrn)
def nullfunc():
	p['na_cond'] = 0
	p['kca_cond'] = 0

 # reduce to one compartment
	#offset = p['iamp']
	# voltage nullcline (cai st. dV/dt = 0
	import nullcline_funcs as nf # this creates new cells
	
	# get holding currents, voltages
	f, varray, iarray,dummy = nf.IV_curve(dopa,args=(0,p),volt_bounds=[-80,10])
	
	ln2 = len(iarray)
	
	#for i in range(ln2):
		#iarray[i]-=offset
	
	fp=open('%svoltage_nullcline.dat' % prefix,'w')
	cell = dopa(0,p)
	area = h.area(0.5,sec=cell.nrn.soma)
	
	na_cond =  550.0e-6	
	kca_cond = 59.0e-6
	
	#find steady state calcium from kca conductance
	for i in range(ln2):
		cond = -iarray[i]/(varray[i]+90) # total conductance
		#if cond < 0:
		#	continue
		#print cond*area
		act = 100*cond/(kca_cond*area)
		#print act
		if act < 0 or act > 1:
			cai=0
		else:
			cai = pow(act,0.25)*190.0/pow(1.0-act,0.25)

		#print varray[i],cai
		fp.write('%e  %e  %e\n' %(varray[i],cai,iarray[i]))
		lastval = varray[i]
	
	fp.write('%e  %e\n' %(lastval,0))
	
	fp.close()
	#quit()
	cells = [cell]
	clamps = []
	va = np.arange(-80,10,1)
	
	i=0
	for voltage in va:
		for s in cells[-1].nrn.all:
			s.fCa_cabalan = 0.02
			s.gnabar_hh3 = 0
		clamps.append(h.SEClamp(cells[-1].nrn.soma(0.5)))
		clamps[-1].dur1 =1e9
		clamps[-1].amp1 = voltage
		clamps[-1].rs=1e-3
		if len(cells)<len(va):
			cells.append(dopa(0,p))
	
	
	cv = h.CVode(1)
	cv.active(1)
	
	h.dt=0.1
	h.finitialize()
	while h.t < 2000: # this is inefficient
		h.fadvance()
		if h.t%100<0.1:
			print h.t
	#cv.solve(10000) # get calcium to steady state
	
	fp = open('%scalcium_nullcline.dat' %prefix,'w')
	for cell in cells:
		fp.write('%e  %e\n' % (cell.nrn.soma.v, 1e6*cell.nrn.soma.cai))
		
		
	fp.close()
	

def runfunc():
	
	cells = []
	cells.append(dopa(0,p))
	cells.append(dopa(0,p))
	
	for s in cells[1].nrn.all:
		s.cm =1e-6
		s.fCa_cabalan=1e-3
		s.gnabar_hh3 = 0
	
	cells.append(dopa(0,p))
	
	for s in cells[2].nrn.all:
		s.gkbar_kca = kca_cond
		s.gnabar_hh3 = 0
	
	
	print cells[2].nrn.soma.gnabar_hh3
	vvecs = []
	cavecs = []
	
	for c in cells:
		vvecs.append(h.Vector())
		cavecs.append(h.Vector())
		vvecs[-1].record(c.nrn.soma(0.5)._ref_v,0.1)
		cavecs[-1].record(c.nrn.soma(0.5)._ref_cai,0.1)
	
	
	
	h.finitialize(p['vinit'])
	h.t = 0
	h.dt = 0.1
	while h.t < 10000:
		h.fadvance()
		if h.t%100 < h.dt:
			print h.t
	print 'done single cell'
	fp = open('%straces.dat' % prefix,'w')
	for i in range(len(vvecs[0])):
		if i < 5000:
			continue
		fp.write('%e  ' % (i*0.1))
		for j in range(3):
			fp.write('%e  %e  ' %(vvecs[j].x[i],1e6*cavecs[j].x[i]))
		fp.write('\n')
			
	fp.close()


cells = None


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
	
	p['kca_cond']=59e-6
	cell = dopa(0,p)

	
	TSTOP=10000
	
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
	
	spikes_gaba = gen_netstim(stepfunc, None,nspikes=10,duration=TSTOP,funcargs={'basal':100,'delay':3000,'duration':1000,'pulseamp':0},ncells=1,syn_type='gabaa')
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
	
	p['kca_cond']=59e-6
	#cell = dopa(0,p)

	nmda_array = np.arange(0,4.01e-5,1e-6)
	gaba_array = np.arange(0,8.01e-5,2e-6)
	ncells = len(nmda_array)*len(gaba_array)
	# 900 cells
	TSTOP=20000


	

	pc=h.ParallelContext()
	pc.nthread(64)
	
	spikes_gaba = gen_netstim(stepfunc, None,nspikes=10,duration=TSTOP,funcargs={'basal':100,'delay':3000,'duration':0,'pulseamp':0},ncells=1,syn_type='gabaa')
	spikes_nmda = gen_netstim(stepfunc, None,nspikes=10,duration=TSTOP,funcargs={'basal':10,'delay':7000,'duration':0,'pulseamp':0},ncells=1,syn_type='nmda')
	cells = []
	gaba_syns = []
	gaba_ncs = []
	nmda_syns= []
	nmda_ncs = []
	i=0
	for nmda_cond in nmda_array:
		for gaba_cond in gaba_array:
			cells.append(dopa(0,p))
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
	fp = open('balarraydata.dat','w')
	for i in range(ncells):
		cells[i].srec.where('>',1000)# chop off initial transient
		nspikes = len(cells[i].srec)
		if nspikes>2:
			isitemp.deriv(cells[i].srec,1)
			freq = 1000.0/isitemp.mean()
			var = isitemp.var()/isitemp.mean()
			burstiness = len(isitemp.where('<',50))/(nspikes-1.0)
		else:
			burstiness = 0
			freq = 0
			var = 0
			
		fp.write('%e  %e  %e  %e  %e\n' %(gaba_ncs[i].weight[0],nmda_ncs[i].weight[0],freq,var,burstiness))
		#print('%e  %e  %e' %(gaba_ncs[i].weight[0],nmda_ncs[i].weight[0],len(cells[i].srec)/(0.001*TSTOP))   )
	fp.close()

if RUN:
	runfunc()
	
if NULLCLINE:
	nullfunc()


if BALANCED:
	balanced()
	
if BALARRAY:
	balanced_array()
quit()
