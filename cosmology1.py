"""
	PROGRAMA PARA CALCULAR DISTANCIAS EM COSMOLOGIA
		Lucas Stefano Xavier de Sousa
			23/12/2020
"""
import numpy as np
import pylab as pl
import scipy.integrate as inte

############## Fiducial values ##################
h=0.72						#
omf =0.25					# Matter
odef = 0.75					# Dark Energy
orf = 8.2E-5					# Radiation
wf = -1.					# Eq. of state
okf = 1. - omf - odef - orf			# Curvature
c = 2.99E5					# c in km/s
H0 = 100*h/c					# Hubble in Km/s/Mpc
#################################################
z1 = np.linspace(0.0,10.,1000)
def H(z,om,ode,orr,ok,w):
	aom = om*np.power((1+z), 3)
	aok = ok*np.power((1+z), 2)
	aor = orr*np.power((1+z), 4)
	aode = ode*np.power((1+z),(3+3*w))
	return 1./(H0*np.sqrt(aom+aok+aor+aode))

def Dx(z,om,ode,orr,ok,w):
	return inte.quad(H,0,z,args=(om,ode,orr,ok,w))[0]
def Dx_matrix(z, om,ode,orr,ok,w):
	Dx_matrixl = np.zeros_like(z)
	for i in range(len(z1)):
		Dx_matrixl[i] = Dx(z[i],om,ode,orr,ok,w)
	return Dx_matrixl


def Da(z,om,ode,orr,ok,w):
	if ok == 0.0:
		return Dx_matrix(z,om,ode,orr,ok,w)
	elif ok < 0.0:
		K = H0*H0*np.abs(ok)
		return np.power((K), -1./2)*np.sin(np.sqrt(K)*Dx_matrix(z,om,ode,orr,ok,w))
	else:
		K = H0*H0*np.abs(ok)
		return np.power((K), -1./2)*np.sinh(np.sqrt(K)*Dx_matrix(z,om,ode,orr,ok,w))

def Dl(z,om,ode,orr,ok,w):
	return Da(z,om,ode,orr,ok,w)*(1+z)



pl.figure()
pl.plot(z1,Dx_matrix(z1,omf,odef,orf,okf,wf))
pl.plot(z1,Da(z1,omf,odef,orf,okf,wf))
pl.plot(z1,Dl(z1,omf,odef,orf,okf,wf))
pl.show()
