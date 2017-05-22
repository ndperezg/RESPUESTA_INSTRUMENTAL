"""
Script para revisar la pruebas de laboratorio de las Q330
solo necesita forma de onda de las pruebas seno
"""
from obspy import read
import sys, os

if len(sys.argv)<2:
	print "No hay parametros suficientes"
	sys.exit()
wave = sys.argv[1]

st = read(wave)
S = 419430.40 # Sensibilidad de Q330 (Count/V)

Vpp = float(raw_input('Ingrese el voltaje pico-pico de la senal de prueba (V):\n'))

def error(V1, V2):
	return abs((V1 - V2)/V1)*100	

for tr in st:
	max_counts = max(tr.data)
	V = max_counts/S
	E = error(Vpp/2.,V)
	if E > 10.0:
		print ("El canal %s tiene problemas: \n Voltaje pico: %4f, Voltaje calculado: %4f, error: %4f %% ")%(tr.stats.channel,Vpp/2., V, E)
	else:
		print ("Canal %s OK: \n Voltaje pico: %4f, Voltaje calculado: %4f, error: %4f %% ")%(tr.stats.channel,Vpp/2., V, E)
		
