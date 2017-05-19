# *-* coding: iso-8859-1 *-*
from obspy.io.xseed import Parser
from obspy import UTCDateTime
from test_dataless import Prod, A0_calc, Transfer_calc, Plot_transfer
import numpy as np
import sys

def rivera(z,p,S1,S2,f):
	Sv = S1*S2
	z.append(complex(0.,0.))      # vel ---> desp

	hzflag=False

	if hzflag:
		for k in xrange(len(z)):
			z[k] *= 2.*np.pi
		for k in xrange(len(p)):
			p[k] *= 2.*np.pi
		
	I = complex(0.,1.)
	R  = 1.
	fr = 1.
	omr = fr * 2 * np.pi

	#R  = C
	for k in xrange(len(z)-1):
		R *= (I*omr - z[k])
	for k in xrange(len(p)):
		R /= (I*omr - p[k])

		
	N = abs(1/R)
	C  = N*Sv

	print 'Factor de Normalización: %15.5e'%N 

	print 'Sensibilidad: %15.5e'%Sv

	#f = arange(.001,100,.001)
	RR = []
	for j in xrange(len(f)):
		om = 2*np.pi*f[j]
		R = C
		for k in xrange(len(z)-1):
			R *= (I*om - z[k])
		for k in xrange(len(p)):
			R /= (I*om - p[k])
		RR.append(R)
	RR, f = np.array(RR), np.array(f)
	return RR, f, Sv, N

PAZ = {'poles': [-0.037004 + 0.037016j, -0.037004 - 0.037016j, -251.33 + 0j,- 131.04 - 467.29j, -131.04 + 467.29j],'zeros': [0j, 0j],'gain': 60077000.0,'sensitivity': 2516778400.0}
F = np.arange(.001,100,.001)
###Calcula con funcion de L. Rivera
RR, f, Sv, N = rivera(PAZ['zeros'], PAZ['poles'],419430.40,1335.00,F)
#Plot_transfer(F,RR,Sv,1.0,1.0,'HHZ','start','end', '00','png')
###Calcula con funciones de test_dataless
H_s = Prod(PAZ['poles'], PAZ['zeros'], 1.0)
N_0 = abs(1/H_s)
A0 = A0_calc(PAZ['poles'], PAZ['zeros'], 1.0)
Sd = 1335.0*419430.40 
RR_, R = Transfer_calc(PAZ['poles'], PAZ['zeros'], F, PAZ['sensitivity'], A0, 1.0)
#Plot_transfer(F, RR_, Sd,1.0, R, 'HHZ_1', 'start', 'end', '00', 'png')
print "sensibilidad %s - %s"%(Sd, Sv)
print "Factor de normalizacion: A0:%s - N:%s - N_0:%s"%(A0, N, N_0)

