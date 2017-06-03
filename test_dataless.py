# *-* coding: iso-8859-1 *-*
from obspy.io.xseed import Parser
from obspy.core.inventory.response import _pitick2latex
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.transforms import blended_transform_factory
import numpy as np
import sys



###Funciones para calcular A0 y H(s) con polos y ceros. 
###Funciona solo para funcion de transferencia tipo A en campo 003 de blockette tipo [53]

def Prod(poles, zeros, fn):
	H_s = 1.
	s = np.complex(0,2.0*np.pi*fn)
	for z in zeros:
		H_s *=(s - z)
	for p in poles:
		H_s /=(s - p)
	return	H_s

def A0_calc(poles, zeros, fn):
	H_s = Prod(poles, zeros, fn)
	return  abs(1/H_s)

def Transfer_calc(poles, zeros, F, Sd, A0, fn):
	RR = []
	C = A0*Sd
	for f in F:
		RR.append(Prod(poles,zeros,f))
		if f == fn:
			R = Prod(poles,zeros,fn)
	return C*np.array(RR), R

def error_perc(A0_d, A0_calc):
	return abs((A0_d - A0_calc)/A0_d)*100

def argand_plot(poles,zeros,channel_id,starttime,endtime):
	poles, zeros = np.array(poles), np.array(zeros)
	fig, ax = plt.subplots()
	ax.scatter(poles.real, poles.imag)
	plt.show()

def Plot_transfer(F, RR, Sd, fn, R, channel_id, starttime, endtime, location, _format):
	lw=1.5
	if location == '10':
		var = 'aceleracion'
		units = 'cuentas/(m/s^2)'
	else:
		var = 'velocidad'
		units = 'cuentas/(m/s)'
		
	xmin = min(F)
	xmax = max(F)
	ymin = min(abs(RR))
	ymax = max(abs(RR)) + 1e1
	xmin,xmax = xmin -(xmax-xmin)/20., xmax+(xmax-xmin)/20.
	ymin,ymax = ymin -(ymax-ymin)/20., ymax+(ymax-ymin)/20.

	fig = plt.figure()

	ax1 = fig.add_subplot(211)
	ax1.loglog(F,abs(RR),lw=lw)
	ax1.scatter(fn,ymax,linewidth=2, s = 30, marker='o',color='b', label='fn = %s Hz'%fn)
	ax1.axvline(x=fn,ymax=max(abs(RR)), linewidth=2, linestyle='--')
	ax1.grid(True)
	ax1.set_xlim([xmin, xmax])
	ax1.set_ylim([ymin, ymax])
	plt.xlabel('')
	plt.ylabel('Sensibilidad %s'%units)
	ax1.xaxis.set_major_formatter(mtick.NullFormatter())
	ax1.yaxis.set_major_locator(mtick.LogLocator(base=10,subs=[1,2,3,4,5,6,7,8,9]))
	ax1.text(xmax - .05*(xmax-xmin), ymin+.05*(ymax-ymin),'Sensibilidad: %13.3e %s'%(Sd,units), size=15,  horizontalalignment='right')
	plt.title('Respuesta en %s: %s %s - %s'%(var,channel_id, starttime, endtime))

	# Phase
	ax2 = fig.add_subplot(212)
	ax2.semilogx(F,np.angle(RR),lw=lw)
	ax2.scatter(fn, 0,linewidth=2, s = 30, marker='o',color='b', label='fn = %s Hz'%fn)
	plt.legend(loc=1,scatterpoints=1, markerscale=1)
	ax2.grid(True)
	ax2.set_xlim([xmin, xmax])
	plt.xlabel('frecuencia [Hz]')
	plt.ylabel('Fase [Radianes]')
	minmax2 = ax2.yaxis.get_data_interval()
   	yticks2 = np.arange(minmax2[0] - minmax2[0] % (np.pi / 2), minmax2[1] - minmax2[1] % (np.pi / 2) + np.pi, np.pi / 2)
    	ax2.set_yticks(yticks2)
    	ax2.set_yticklabels([_pitick2latex(x) for x in yticks2])
	name_fig = '%s_%s.%s'%(channel_id,starttime,_format)
	plt.savefig(name_fig,dpi=fig.dpi)
	print "%s saved"%name_fig
	plt.show()
	
if __name__ == '__main__':
	#formato de las imagenes
	_format = 'png'

	#Parametros de entrada			
	dataless_name = sys.argv[1]
	#Frecuencia de normalizacion. Para estados analogos se toma normalmente como  1.0 Hz
	fn = float(raw_input("Ingrese el valor de la frecuencia de referencia en Hz\n (Valor recomendado: 1.0 Hz)\n")) 

	#Crea archivo log
	log_file = open(dataless_name+'.log', 'w')

	#Carga dataless mediante io.xseed
	dtlss = Parser(dataless_name)
	print >> log_file, dtlss

	#Crea un diccionario con respuesta por canal
	inv = dtlss.get_inventory()

	###workout para el dataless:
	F = np.arange(.001,100,.001)
	for chn in inv['channels']:
		channel_id, start_date, end_date, instrument = chn['channel_id'], chn['start_date'], chn['end_date'], chn['instrument']
		location = channel_id.split('.')[2]
		if end_date != "":
			starttime, endtime = start_date.date, end_date.date
		else:
			starttime, endtime = start_date.date, " "
		PAZ = dtlss.get_paz(seed_id=channel_id,datetime=start_date)
		A0 = A0_calc(PAZ['poles'], PAZ['zeros'], fn)
		Sd = PAZ['seismometer_gain']*PAZ['digitizer_gain'] 
		RR, R = Transfer_calc(PAZ['poles'], PAZ['zeros'], F, PAZ['sensitivity'], A0, fn)
		#argand_plot(PAZ['poles'], PAZ['zeros'],channel_id,starttime,endtime)
		Plot_transfer(F, RR, Sd,fn, R, channel_id, starttime, endtime, location, _format)
		print >> log_file, 100*'='+'\n'
		print >> log_file, channel_id+'  '+str(start_date)+'  '+str(end_date)+'\n'
		print >> log_file, "A0 dataless: %15.5e     A0 calculado: %15.5e   error(%%): %6f \n"%(PAZ['gain'], A0, error_perc(PAZ['gain'],A0)) 	
		print >> log_file, "Sensibilidad dataless: %15.5e    Sensibilidad calculada: %15.5e  error(%%): %6f\n"%(PAZ['sensitivity'], Sd, error_perc(PAZ['sensitivity'], Sd))
	###Imprime en pantalla:

	log_file.close()
	for line in open(dataless_name+'.log','r').readlines():
		sys.stdout.write(line)
