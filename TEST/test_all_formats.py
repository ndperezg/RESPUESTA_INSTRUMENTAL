# -*- coding: utf-8 -*-
"""
Programa que revisa archivos de respuesta en formato Dataless, GSE, SAC e ISOLA
v(0.0)
autor Nelson David Perez - nperez@sgc.gov.co
"""

import argparse
import sys, os
import numpy as np
from obspy.io.xseed import Parser
from obspy.core.inventory.response import _pitick2latex
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from test_dataless import A0_calc, Transfer_calc


def read_GSE(GSE):
	with open(GSE) as gse:
		lines_org, poles, zeros = [], [], []
		for line in gse:
			split = line.split()
			if len(split)==2:
				lines_org.append(line.strip())
			if 'CAL2' in split:
				Gtot = float(split[4])
				label = ".".join((split[1],split[2]))
			if 'DIG2' in split:
				Gd = float(split[2])
			if 'PAZ2' in split:
				npoles, nzeros = int(split[4]), int(split[5])
		for k in xrange(npoles):
 		         poles.append(complex(float(lines_org[k].split()[0]), float(lines_org[k].split()[1])))
		for k in xrange(nzeros):
			zeros.append(complex(float(lines_org[k+npoles].split()[0]), float(lines_org[k+npoles].split()[1])))
		Gs = (1e9)/(2*np.pi*Gtot*Gd)
		print label	
	return poles, zeros, Gs*Gd, label
			
def read_SEED(SEED):
	dataless = Parser(SEED)
	inv = dataless.get_inventory()
	chn = list(inv['channels'])
	for i in range(len(chn)):
		channel_id, start_date, end_date, instrument = chn[i]['channel_id'], chn[i]['start_date'], chn[i]['end_date'], chn[i]['instrument']
		print "Channel %s: \t %s \t %s \t %s"%(i, channel_id, start_date, end_date)
	n = int(raw_input("Ingrese el número del canal que quiere probar\n"))
	channel_id, start_date, end_date, instrument = chn[n]['channel_id'], chn[n]['start_date'], chn[n]['end_date'], chn[n]['instrument']
	PAZ = dataless.get_paz(seed_id=channel_id,datetime=start_date)
	label = ".".join((channel_id.split('.')[1],channel_id.split('.')[3]))
	print label
	poles, zeros, Sd = PAZ['poles'], PAZ['zeros'], PAZ['sensitivity']
	return poles, zeros, Sd, label
		 		

def read_SAC(SAC):
	with open(SAC) as sac:
		label = ".".join((sac.name.split('_')[3],sac.name.split('_')[4]))
		print label
		Sd, poles, zeros = [], [], []
		lines = sac.readlines()
		for i in xrange(len(lines)):
			split = lines[i].split()
			if 'ZEROS' in split:
				nz = int(split[1])
				kz = i+1
			if 'POLES' in split:
				np = int(split[1])
				kp = i+1
			if 'CONSTANT' in split:
				C = float(split[1])
		for j in lines[kz+1:kz+nz]:
			zeros.append(complex(float(j.strip().split()[0]),float(j.strip().split()[1])))
		for l in lines[kp:kp+np]:
			poles.append(complex(float(l.strip().split()[0]),float(l.strip().split()[1])))
		A0 = A0_calc(poles,zeros,fn=1.0)
		return poles, zeros, C/A0, label
 
def read_Isola(ISOLA):
	with open(ISOLA) as isola:
		label = ".".join((isola.name.split('.')[0][:-3],isola.name.split('.')[0][-3:]))
		print label
		org = isola.readlines()
		numbers, poles, zeros = [], [], []
		for i in xrange(len(org)):
			if org[i].strip() == 'count-->m/sec':
				Sd = 1./float(org[i+1])
			if org[i].strip() == 'zeroes':
				kz = int(org[i+1])
			if org[i].strip() == 'poles':
				kp = int(org[i+1])
			if len(org[i].strip().split()) == 2:
				numbers.append(org[i].strip())
		for j in xrange(kz):
			zeros.append(complex(float(numbers[j].split()[0]),float(numbers[j].split()[1])))
		for l in xrange(kp):
			poles.append(complex(float(numbers[kz+l].split()[0]),float(numbers[kz+l].split()[1])))
		return poles, zeros, Sd, label
			

def dirty(File, f):
	if f == 'gse':
		try:
			marker = 'o'
			poles, zeros, Sd, label = read_GSE(File)
			print "The file %s was read"%File
		except:
			print "error in file %s"%File
	if f == 'sac':
		try:
			marker = '^'
			poles, zeros, Sd, label = read_SAC(File)
			print "The file %s was read"%File
		except:
			print "error in file %s"%File
	if f == 'isola':
		try:
			marker = 'd'
			poles, zeros, Sd, label = read_Isola(File)
			print "The file %s was read"%File
		except:
			print "error in file %s"%File
	if f == 'dataless':
		try:
			marker = 'o'
			poles, zeros, Sd, label = read_SEED(File)
			print "The file %s was read"%File
		except:
			print "error in file %s"%File
	return poles, zeros, Sd, label, marker 

def bode_plot_generator(F):
	xmin = min(F)
	xmax = max(F)
	ymin = 1e7
	ymax = 1e10 
	fig = plt.figure(figsize=(10,10))
	ax1 = fig.add_subplot(211)
	ax1.grid(True)
	plt.xlabel('')
	plt.ylabel('Sensibilidad')
	#ax1.set_xlim([xmin, xmax])
	#ax1.set_ylim([ymin, ymax])
	ax1.xaxis.set_major_formatter(mtick.NullFormatter())
	ax1.yaxis.set_major_locator(mtick.LogLocator(base=10,subs=[1,2,3,4,5,6,7,8,9]))

	# Phase
	ax2 = fig.add_subplot(212)
	ax2.grid(True)
	plt.xlabel('frecuencia [Hz]')
	plt.ylabel('Phase [Radianes]')
	return fig, ax1, ax2

def pi_ticks(ax2):
	minmax2 = ax2.yaxis.get_data_interval()
   	yticks2 = np.arange(minmax2[0] - minmax2[0] % (np.pi / 2), minmax2[1] - minmax2[1] % (np.pi / 2) + np.pi, np.pi / 2)
    	ax2.set_yticks(yticks2)
    	ax2.set_yticklabels([_pitick2latex(x) for x in yticks2])

def arg_generator(parser,Dict):
	for key in Dict.keys():
		parser.add_argument(key, dest=Dict[key], metavar="FILE", nargs="+", help="file in %s format"%(Dict[key]))
	args = parser.parse_args()
	return args

def main():
	choices={"-g":"gse", "-d":"dataless", "-s":"sac", "-i":"isola"}
	parser = argparse.ArgumentParser(description='This program plots Amplitude and frequency response for Dataless SEED, GSE, SAC or ISOLA response files')
	args = arg_generator(parser, choices)
	files = vars(args)
	F = np.arange(.001,100,.001)
	fn = 1.
	fig, ax1, ax2 = bode_plot_generator(F)
	lw=1.5
	mew=0.15
	counter = 0
	for f in files.keys():
		lst = files[f]
		if lst != None:
			counter += 1 
			for i in lst:
				#print i, f
				poles, zeros, Sd, label, marker = dirty(i,f)
				A0 = A0_calc(poles,zeros,fn)
				RR, R = Transfer_calc(poles, zeros, F, Sd, A0, fn)
				ax1.loglog(F,abs(RR),lw=lw, marker=marker, mew=mew,alpha = 0.3)
				ax2.semilogx(F,np.angle(RR),lw=lw, marker=marker,mew=mew, alpha = 0.3)
		else: 
			pass
	if counter > 0:
		pi_ticks(ax2)	 
		plt.show()

if __name__ == "__main__":
	print "Running %s"%(__file__)
	main()
