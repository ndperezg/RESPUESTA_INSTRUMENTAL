# -*- coding: utf-8 -*-
"""
Programa que revisa archivos de respuesta en formato Dataless, GSE, SAC e ISOLA
v(0.0)
autor Nelson David Péez - nperez@sgc.gov.co
"""

import argparse
import sys, os

def read_GSE(GSE):
	with open(GSE) as gse:
		for line in gse:
			print line

choices=['g','d','s','i']

def main():
	parser = argparse.ArgumentParser(description='This program plots Amplitude and frequency response for Dataless SEED, GSE, SAC or ISOLA response files')
	parser.add_argument("-g", dest='gse', help='File in GSE format',action="store_true")
	parser.add_argument("-d", dest='dataless', help='File in Dataless SEED format',action="store_true")
	parser.add_argument("-s", dest='sac', help='File in Dataless SEED format',action="store_true")
	parser.add_argument("-i", dest='isola', help='File in Dataless SEED format',action="store_true")
	parser.add_argument("filename", help= "Response file name")
	args = parser.parse_args()
	if args.gse:
		read_GSE(args.filename)
	else:
		print 'OK' 

if __name__ == "__main__":
	print "Running %s"%(__file__)
	main()
