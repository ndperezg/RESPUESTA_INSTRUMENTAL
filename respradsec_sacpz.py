#!/usr/bin/env python
# *-* coding: iso-8859-1 *-*

import argparse
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from sys import stdout, argv

parser = argparse.ArgumentParser()
parser.add_argument("Input_File", 
help='''Archivo de entrada de respuesta en hz.  
        Este script prepara un archivo de respuesta
        en polos y ceros en formato sac a partir de 
        un archivo de entrada con informacion de 
        respuesta en rad/s con formato ad hoc''')
args = parser.parse_args()

print args.Input_File

# Input
ifile = argv[1]
lines_org = open(ifile,'r').readlines()
lines = []
for k in xrange(len(lines_org)):
	if lines_org[k][0] != '#':
		lines.append(lines_org[k].split('#')[0])

# Decoding
net,stat,chan,loc = lines[0].split()
S2     = float(lines[1])
S1     = float(lines[2])
npoles =   int(lines[3])
p = []
for k in xrange(npoles):
	p.append(complex(float(lines[4+k].split()[0]), float(lines[4+k].split()[1])))

#print  npoles, p

k0 = k+5
nzeros = int(lines[k0])
print nzeros
z = []
for k in xrange(nzeros):
	z.append(complex(float(lines[k0+k+1].split()[0]), float(lines[k0+k+1].split()[1])))
print nzeros, z

##############################################
Sv = S1*S2
z.append(complex(0.,0.))      # vel ---> desp

hzflag=False

if hzflag:
	for k in xrange(len(z)):
		z[k] *= 2.*pi
	for k in xrange(len(p)):
		p[k] *= 2.*pi
	
I = complex(0.,1.)
C  = 1.
fr = 1.
omr = fr * 2 * pi

R  = C
for k in xrange(len(z)-1):
	R *= (I*omr - z[k])
for k in xrange(len(p)):
	R /= (I*omr - p[k])

	
N = abs(1/R)
C  = N*Sv

print 'Factor de Normalización: %15.5e'%N

print 'Sensibilidad: %15.5e'%Sv

f = arange(.001,100,.001)
RR = []
for j in xrange(len(f)):
	om = f[j]
	R = C
	for k in xrange(len(z)-1):
		R *= (I*om - z[k])
	for k in xrange(len(p)):
		R /= (I*om - p[k])
	RR.append(R)
#	print abs(R)

fname = '%s_%s_%s_%s'%(net,stat,chan,loc)
sac_pz = 'SAC_PZs_'+fname
pdf    = 'resp_'+fname+'.pdf'
fd = open(sac_pz,'w')

fd.write('ZEROS %d\n'%len(z))
for k in xrange(len(z)):
	fd.write('%15.5e  %15.5e\n'%(real(z[k]), imag(z[k])))
fd.write('POLES %d\n'%len(p))
for k in xrange(len(p)):
	fd.write('%15.5e  %15.5e\n'%(real(p[k]), imag(p[k])))

fd.write('CONSTANT %e\n'%C)

fd.close()

# Console output
for line in open(sac_pz,'r').readlines():
	stdout.write(line)

print fname, pdf

# Graphics
f  = array(f)
RR = array(RR)

xmin = min(f)
xmax = max(f)
ymin = min(abs(RR))
ymax = max(abs(RR))
xmin,xmax = xmin -(xmax-xmin)/20., xmax+(xmax-xmin)/20.
ymin,ymax = ymin -(ymax-ymin)/20., ymax+(ymax-ymin)/20.

fig = plt.figure()

ax1 = fig.add_subplot(211)
ax1.loglog(f,abs(RR))
ax1.grid(True)
ax1.set_xlim([xmin, xmax])
ax1.set_ylim([ymin, ymax])
plt.xlabel('')
plt.ylabel('Sensibilidad [cuentas/(m/s)]')
ax1.xaxis.set_major_formatter(mtick.NullFormatter())
ax1.yaxis.set_major_locator(mtick.LogLocator(base=10,subs=[1,2,3,4,5,6,7,8,9]))
ax1.text(xmax - .05*(xmax-xmin), ymin+.05*(ymax-ymin),'Sensibilidad: %13.3e cuentas/(m/s)'%Sv, size=15,  horizontalalignment='right')
plt.title('Respuesta en velocidad : %s %s %s %s'%(net, stat, chan, loc))

# Phase
ax2 = fig.add_subplot(212)
ax2.semilogx(f,angle(RR))
ax2.grid(True)
ax2.set_xlim([xmin, xmax])
ax2.set_ylim([-pi, pi])
ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%5.1f'))
plt.xlabel('frecuencia [Hz]')
plt.ylabel('Phase [Radianes]')

plt.savefig(pdf)
plt.show()
