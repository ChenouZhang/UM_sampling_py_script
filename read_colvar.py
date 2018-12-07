import numpy as np
import matplotlib.pyplot as plt
import math

#Python scripts to read,plot 2D_umbrella_sampling result.
#Add WHAM-2d read and plot.

def ReadColvar_Coord(f):
	'''
	Parameters: f: Input file head
	               Read the x,y coordinates directly from plumed ouput dat file.
	               #! FIELDS time p1.sss p1.zzz restraint-phi.bias restraint-psi.bias
	                         0.00   1.11  0.001       1.105              0.0002

	Returns:    Coord:  array_like coordinates (x,y) pairs
	            x: array_like x coordinates
	            y: array_like y coordinates
	'''
	print(f)
	## Skip first row.
	next(f)
	coord = []
	x = []
	y = []
	for line in f:
		temp1 = float(line.split()[1])
		temp2 = float(line.split()[2])
		x.append(temp1)
		y.append(temp2)
		coord.append([temp1,temp2])
		#y.append(float(temp))
	return coord,x,y

def ReadColvar_Bias(f):
	'''
	Parameters: 
		f: Input file head
			Read the bias potential directly from plumed ouput dat file.
			#! FIELDS time p1.sss p1.zzz restraint-phi.bias restraint-psi.bias
	                  0.00   1.11  0.001       1.105              0.0002

	Returns:
		Bias:  array_like bias_energy (phi_bias,psi_bias) pairs
		phi_bias: array_like x coordinates
		psi_bias: array_like y coordinates
	'''
	print(f)
	## Skip first row.
	next(f)
	Bias = []
	phi_bias = []
	psi_bias = []
	for line in f:
		temp3 = float(line.split()[3])
		temp4 = float(line.split()[4])
		phi_bias.append(temp3)
		psi_bias.append(temp4)
		Bias.append([temp3,temp4])
		#y.append(float(temp))
	return Bias,phi_bias,psi_bias

def ReadWHAM(f):
	'''
	Parameters:
	f: Input file head
		Read free energy and coordinates from WHAM output
		#X      Y      Free     Pro
		0.00   1.11   0.001    0.00

	Returns:
		Coord: array_like coordinates (x,y,Free) pairs
		Free: array_like Free energy
		x: array_like x coordinates
		y: array_like y coordinates
	'''
	print(f)
	## Skip first row
	next(f)
	coord = []
	#x = []
	#y = []
	Free = []
	## Line counter:  WHAM output skip 1 line after 50 lines
	Line_counter = 0
	for line in f:
		if Line_counter == 50:
			Line_counter = 0
			continue
		temp1 = float(line.split()[0])
		temp2 = float(line.split()[1])
		temp3 = float(line.split()[2])
		#x.append(temp1)
		#y.append(temp2)
		Free.append(temp3)
		coord.append([temp1,temp2,temp3])
		#print ('ok'+str(i))
		Line_counter = Line_counter + 1
	return coord,Free

def PlotWHAM_2D(f,binx=1000,biny=50,minx=1.,maxx=11.,miny=-0.005,maxy=0.06):
	'''
	Read WHAM-2D output and generate 'contour' ready outputs

	Parameters: f: Input file head
	            binx: 
				biny:
				minx:
				maxx:
				miny:
				maxy:
	Returns:    H:
				xv:
				yv:
	'''
	x = []
	y = []
	## Set x,y interval
	xint = (maxx-minx)/binx
	yint = (maxy-miny)/biny
	## Get 2D meshgrid on x,y
	x = np.linspace(minx,maxx,binx)
	y = np.linspace(miny,maxy,biny)
	xv,yv = np.meshgrid(x,y)
	## 
	H = [[0]*binx for _ in range(biny)]
	print ('Following are the input parameters:_______________________________________________________________')
	print ('X,Y interval have been set to:'+str(xint)+' '+str(yint))
	print ('binx:'+str(binx),'biny:'+str(biny))
	#print (np.shape(H))
	#xtemp = []
	#ytemp = []
	Coord_temp,Free = ReadWHAM(f)
	num = np.shape(Coord_temp)
	print (num)
	for k in Coord_temp:
		xind = int(math.ceil((k[0]-minx)/xint)-1)
		yind = int(math.ceil((k[1]-miny)/yint)-1)
		if k[2] >= 1000000:
			H[yind][xind] = np.NAN
		else:
			H[yind][xind] = k[2]				
	return H,xv,yv

def plotcolvar(x,y,binx,biny):
	## 2D- histogram plot (Use contourcolvar instead)
	fig = plt.figure()
	plt.hist2d(x,y,bins=(binx,biny))
	plt.xlabel('s')
	plt.xlabel('z')
	cbar = plt.colorbar()
	cbar.ax.set_ylabel('Counts')
	plt.show()
	return 0

def contourcolvar(x,y,binx,biny):
	H,xedge,yedge = np.histogram2d(x,y,bins=(binx,biny))
	xedge = 0.5*(xedge[:-1] + xedge[1:])
	yedge = 0.5*(yedge[:-1] + yedge[1:])
	#CS = plt.contourf(yedge,xedge,H,30,cmap=plt.cm.bone)
	#plt.show()
	#print xdim,ydim,Hdim
	return xedge,yedge,H

def ReadWindows(window_sn,window_zn,binx,biny,minx=1.,maxx=11.,miny=-0.005,maxy=0.07):
	## Input file name as index_sn(1-11) , index_zn(1,2,3)
	## Bin number binx and biny
	##_________________________________________________________________________________________________________
	## Return H,x,y ready for any contour plot.
	x = []
	y = []
	print ('Following are the input parameters:_______________________________________________________________')
	print ('window_sn:'+str(window_sn),'window_zn:'+str(window_zn),'binx:'+str(binx),'biny:'+str(biny))
	xint = (maxx-minx)/binx
	yint = (maxy-miny)/biny
	print ('X,Y interval have been set to:'+str(xint)+' '+str(yint))
	x = np.linspace(minx,maxx,binx)
	y = np.linspace(miny,maxy,biny)
	xv,yv = np.meshgrid(x,y)
	H = [[0]*binx for _ in range(biny)]
	#print (np.shape(H))
	#xtemp = []
	#ytemp = []
	for i in window_sn:
		if i == maxx:
			continue
		for j in window_zn:
			## -----------------------------------------------------------
			if i == 2 and j == 3:
				continue
			## -----------------------------------------------------------
			with open('colvar_reference_'+str(i)+'_'+str(j)+'_20ns') as f:
				Coord_temp,xtemp,ytemp = ReadColvar_Coord(f)
				num = np.shape(Coord_temp)
				print (num)
				for k in Coord_temp:
					xind = int(math.ceil((k[0]-minx)/xint)-1)
					yind = int(math.ceil((k[1]-miny)/yint)-1)
					if k[1]<miny or k[1]>maxy:
						print (k[1])
						continue
					elif k[0]<minx or k[0]>maxx:
						continue
					H[yind][xind] = H[yind][xind]+1
					#if xind >= 1000 or xind <= 0:
					#	print(k)
					#if yind >= 50 or yind <= 0:
					#	print(k)
					#print(xind,yind)
				#print (num)
				
	return H,xv,yv

def FindFWHM(histogram):
	difference = max(histogram) - min(histogram)
	H = np.array(histogram)
	nearest = np.abs(H - difference/2)
	fmin = nearest[0]
	smin = nearest[0]
	fminn = 0
	sminn = 0
	n = 0
	for i in nearest:
		n = n+1
		if i < fmin:
			smin = fmin
			fmin = i
			sminn = fminn
			fminn = n
		elif i > fmin and i < smin:
			smin = i
			sminn = n
	#print (nearest)
	FWHM = np.abs(fminn-sminn)
	return FWHM

##