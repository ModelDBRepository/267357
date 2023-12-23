# functions for creating bifurcation diagrams and nullcline pictures
# of NEURON models

import numpy as np
from neuron import h
import neuron_jacobian as nj
import gc
from mpi4py import MPI # for h.ParallelContext
from scipy import optimize, interpolate

VERBOSE = 0

def IV_curve(nrn_class, args=None, cells=[], volt_bounds=[-80,-20], voltstep=0.25,rundur=10000): #nrn class assumed to have form in dep_fun.py
	
	v = min(volt_bounds)
	stimvec = []
	i=0

	h.t = 0
	h.dt = 0.2
	pc = h.ParallelContext(64)
	pc.nthread(64)
	keepcells=1
	
	if len(cells)>0:
		if not isinstance(cells[0],nrn_class):
			cells = []
		else:
			keepcells = 0
			
	while v <= max(volt_bounds):

		if len(cells) < i+1:
			#h.pop_section()
			#if i > 50:
			#	break
			#h.topology()
			cells.append(nrn_class(*args))
		#for s in cells[i].nrn.:
		stimvec.append(h.SEClamp(cells[i].nrn.soma(0.5)))
		stimvec[i].amp1 = v
		stimvec[i].dur1 =1e9
		stimvec[i].rs = 1e-6
		#h.pop_section()
		#print cells[i].nrn.soma.dist_NaMark
		#varray.append(cells[i].nrn.soma.v)
		v+=voltstep
		i+=1
	
	h.finitialize()
	pc.psolve(rundur)
	varray = []
	iarray = []
	i=0
	for stim in stimvec:
		#print stim.i, cells[i].nrn.soma.v
		iarray.append(stim.i)
		#for s in cells[i].nrn.all:
		varray.append(cells[i].nrn.soma(0.5).v)
		#	break
		i+=1
	
	f = interpolate.InterpolatedUnivariateSpline(varray,iarray,k=3)
	#f=0
	if keepcells:
		return f, varray, iarray, cells
	else:
		return f, varray, iarray, []
		
def VI_curves(nrn_class, args=None, cells=[], volt_bounds=[-80,-20], voltstep=0.25):
	#print args
	f, varray, iarray, cells = IV_curve(nrn_class, args=args, cells=cells, volt_bounds=volt_bounds, voltstep=voltstep)
	fp = f.derivative(1)
	fpp = f.derivative(2)
	step = 10
	v = min(volt_bounds)
	extrema = []
	fpl = fp(varray)
	fppl = fpp(varray)
	doit = True
	#print min(fpl), max(fpl), min(fppl), max(fppl), min(varray), max(varray)
	if not min(fpl) < 0 < max(fpl):
		doit = False
	

	
	while v < max(volt_bounds) and doit:
		temp = optimize.fsolve(fp,v)
		if len(extrema) == 0 and len(temp) > 0:
			extrema.append(temp[0])
		v += step

		for things in temp:
			isin=False
			for others in extrema:
				if abs(things-others) < 1e-1:
					isin = True
			if not isin:
				extrema.append(things)
	extrema.sort()		
	

	extrema = [things for things in extrema if min(volt_bounds) < things < max(volt_bounds)]
	iextrema = f(extrema)
	finverses = []			
	if len(extrema) == 0:
		finverses.append(interpolate.InterpolatedUnivariateSpline(iarray,varray,k=3))
	else:
		extremeold = min(volt_bounds)
		for extreme in extrema:
			#print extreme
			vtemp = [things for things in varray if extremeold <= things <= extreme]
			vtemp.append(extreme)
			itemp = f(vtemp)
			if itemp[-1] < itemp[0]:
				#itemp.reverse()
				temp = itemp[::-1]
				itemp = temp
				temp = vtemp[::-1]
				vtemp = temp
			try:
				finverses.append(interpolate.InterpolatedUnivariateSpline(itemp,vtemp,k=3))
			except:
				pass
			extremeold=extreme
		
	
		vtemp = [things for things in varray if extremeold <= things <= max(volt_bounds)]
		vtemp.append(max(volt_bounds))
		itemp = f(vtemp)
		if itemp[-1] < itemp[0]:
			temp = itemp[::-1]
			itemp = temp
			temp = vtemp[::-1]
			vtemp = temp
		finverses.append(interpolate.InterpolatedUnivariateSpline(itemp,vtemp,k=3))
	
	print finverses, extrema, iextrema
	return finverses, extrema, iextrema, cells


def v_root(nrn_class, cells=[],args=None, init=-40, volt_bounds=[-80,-30],voltstep=1):
	if len(cells) > 0:
		keepcells = 1
	f, varray, iarray, cells = IV_curve(nrn_class,volt_bounds=volt_bounds,args=args,cells=cells,voltstep=voltstep)
	nroots = 0
	
	
	vinits = []
	for i in range(len(iarray)-1):
		if iarray[i]*iarray[i+1] <0:
			nroots+=1
			vinits.append(varray[i])
	#print nroots
	zeros = []
	if nroots > 0:
		for v in vinits:
			fp = f.derivative(1)
			zero = optimize.fsolve(f,v,fprime=fp)
			for z in zero:
				zeros.append(z)
	
	return zeros, cells


def clear_cells(nrn_class, nuke=False):
	#for sec in h.allsec():
		#h.delete_section(sec=sec) # this kills ALL the neuron objects
	if nuke: # this kills all things in the class, probably slow.
		for obj in gc.get_objects():
			if isinstance(obj, nrn_class):
				for s in obj.nrn.all:
					h.delete_section(sec=s)
				del obj

def stability_check(nrn_class, voltage, args=None,offset=0, rundur=1000):
	# need to ensure there are no objects created
	h.t = 0
	h.dt = 0.2
	clear_cells(nrn_class,nuke=True) # clear cells may be required

	cell = nrn_class(*args)

	#print cell.nrn.soma.cai
	cv = h.CVode()
	cv.active(1)
	h.t = 0
	h.dt = 0.1
	if offset != 0: # for stability under iclamp, though it should not affect jacobian
		stim = h.IClamp(cell.nrn.soma(0.5))
		stim.dur = 1e9
		stim.amp = offset
		stim.delay = 0
		
	vc = h.SEClamp(cell.nrn.soma(0.5))
	vc.amp1 = voltage
	vc.dur1 = 1e9
	vc.rs = 1e-6
	
	stat = h.SaveState()
	h.finitialize(voltage)
	cv.solve(rundur)
	#print cell.nrn.soma.cai
	stat.save()
	vc = None # need to kill the voltage clamp as it 'counts' during the jacobian calculation
			  # this could actually be very useful for find when voltage clamp is ineffective
			  # see literature on control theory, coupled oscillators
	
	# This is a kludge - ideally there would be a check, then warning indicating
	# that one or more states did not get copied via savestate/restore
	# 
	for s in cell.nrn.all:  
		if h.ismembrane('cabalan',sec=s):
			s.cainit_cabalan = s.cai
	
	#nj.get_size(cv)	
	h.finitialize(voltage)

	
	stat.restore(1)
	#h.fcurrent()
	#nj.get_size(cv)
	#print cell.nrn.soma.v
	
	#sz = nj.get_size(cv)
	#print sz
	#quit()
	jac_obj = nj.get_jacobian(cv,relstep=1e-6)
	"""if voltage > -35:
		sz = nj.get_size(cv)
		for lines in jac_obj:
			for items in lines:
				print items,
			print
		#print cell.nrn.soma.gkbar_kca
		quit()#"""
	eigs, eigvec = np.linalg.eig(jac_obj)
	
	#print eigs
	
	max_eig = -1e9
	ncomp = 0
	
	for i in range(len(eigs)):
		#print np.imag(eigs[i]),
		if abs(np.imag(eigs[i])) < 1e-12:
			eigs[i] = np.real(eigs[i])
			if abs(eigs[i]) < 1e-12:
				eigs[i] = 0 
		else:
			ncomp +=1
	#print
	comp = 0
	for things in eigs:
		if np.real(things) > max_eig and abs(np.real(things)) > 1e-12:
			max_eig = np.real(things)
			comp = np.imag(things)
	del cv
	del cell
	return max_eig, comp, ncomp, eigs
