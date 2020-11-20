import numpy as np
from math import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def graph(XE,YE,XRK2,YRK2,XRK4,YRK4,title,xlab='x',ylab='y',clf=True,log=False):
	if clf:
		plt.clf()
	fig, ax = plt.subplots()
	#plt.plot(XX,YX, '-',linewidth=3)
	ax.plot(XE,YE, 'o-',linewidth=3)
	ax.plot(XRK2,YRK2, '--',linewidth=3)
	ax.plot(XRK4,YRK4, ':',linewidth=3)
	plt.title(title)#+" "+ylab+"("+xlab+")")
	plt.xlabel(xlab); plt.ylabel(ylab); ylab = ylab +'\''
	xmin = min(min(XE),min(XRK2),min(XRK2))
	xmax = max(max(XE),max(XRK2),max(XRK2))
	plt.xlim(xmin,xmax)
	plt.legend(["Euler","RK2","RK4"])
	if log:
		ax.set_yscale("log")
		#locmaj = ticker.LogLocator(base=10,numticks=15) 
		#ax.yaxis.set_major_locator(locmaj)
		#locmin = ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=15)
		#ax.yaxis.set_minor_locator(locmin)
		#ax.yaxis.set_minor_formatter(ticker.NullFormatter())

	plt.grid(True, which='both')
	plt.gcf().set_size_inches(18.5, 10.5)
	plt.savefig(title+'.png',bbox_inches='tight')


def graphX(XX,YX,X,Y,title,xlab='x',ylab='y',clf=True):
	if clf:
		plt.clf()
	for i in range(2):
		plt.subplot(1,2,i+1)
		plt.plot([X,Y][i],[Y,X][i],'o-',linewidth=3)
		plt.plot([XX,YX][i],[YX,XX][i],linewidth=3)
		plt.title(title+" "+ylab+"("+xlab+")")
		plt.grid(True)
		plt.xlabel(xlab); plt.ylabel(ylab);
		plt.xlim(min(min([XX,YX][i]),min([X,Y][i])),max(max([XX,YX][i]),max([X,Y][i])))
		plt.legend([title,"Exact"])
		plt.gcf().set_size_inches(18.5, 10.5)
		plt.savefig(title+'.png',bbox_inches='tight')
		xlab,ylab = ylab,xlab


#Question 1
#Part A
#y_{i+1}=y_i+h*y'_i=y_i+h*F(x,y)
def euler(F, x, y, xStop, h):
	X = np.array(x)
	Y = np.array(y)
	while x < xStop:
	  #Avoid overstepping
		h = min(h, xStop - x)
		#y'=F(x,y)
		y = y + h*F(x,y)#Use equation 6
		#to calculate the next y
		x = x + h
		X = np.append(X, x)
		Y = np.append(Y, y)
	return X,Y

#y(x+h)=y(x)+F[x+h/2+y+h/2*F(x,y)]h
def RK2(F, x, y, xStop, h):
	X = np.array(x)
	Y = np.array(y)
	while x < xStop:
		h = min(h, xStop - x)#dont overstep
		k0 = h*F(x,y)#eq 14
		k1 = h*F(x + h/2, y + k0/2)#eq 15
		y = y + k1# eq 16
		x = x + h#increment x
		X = np.append(X, x)
		Y = np.append(Y, y)
	return X,Y

def RK4(F, x, y, xStop, h):
	X = np.array(x)
	Y = np.array(y)
	while x < xStop:
		h = min(h, xStop - x)#dont overstep
		k0 = h*F(x,y)#eq 20
		k1 = h*F(x + h/2, y + k0/2)#eq 21
		k2 = h*F(x + h/2, y + k1/2)#eq 22
		k3 = h*F(x + h, y + k2)#eq 23
		y = y + (k0+2*k1+2*k2+k3)/6#eq 24
		x = x + h#increment x for next iteration
		X = np.append(X, x)
		Y = np.append(Y, y)
	return X,Y

def Q2A():
	def F(x,y):
		return sin(y)
	x, xStop = 0.0, 2.0
	y, h = 1.0, 0.1
	XE,YE = euler(F,x,y,xStop,h)
	X2,Y2 = RK2(F,x,y,xStop,h)
	X4,Y4 = RK4(F,x,y,xStop,h)
	graph(XE, YE, X2, Y2, X4, Y4, "y'=sin(y(x))")
	return XE,YE,Y2,Y4

def Q2B(XE,YE,Y2,Y4):
	def X(x,y):
		return np.log(1/np.sin(y)-1/np.tan(y)) + 0.604582
	graphX(YE,X(YE,YE),YE,XE,"Euler",'y','x')
	graphX(Y2,X(Y2,Y2),Y2,XE,"Runge-Kutta 2",'y','x')
	graphX(Y4,X(Y4,Y4),Y4,XE,"Runge-Kutta 4",'y','x')
	return

def Q2C(XE,YE,Y2,Y4):
	def X(x,y):
		return np.log(1/np.sin(y)-1/np.tan(y)) + 0.604582
	XXE = X(YE,YE)	
	XX2 = X(Y2,Y2)	
	XX4 = X(Y4,Y4)
	XEerr = XE - XXE
	X2err = XE - XX2
	X4err = XE - XX4
	graph(XE,XEerr,XE,X2err,XE,X4err,"Errors",'y','x(y) err')
	#graph(YE,XEerr,Y2,X2err,Y4,X4err,"Errors",'y','x(y) err')
	
	XEerr2 = XEerr**2
	X2err2 = X2err**2
	X4err2 = X4err**2
	print("Euler: ",np.sum(XEerr2))
	print("RK2: ",np.sum(X2err2))
	print("Rk4: ",np.sum(X4err2))
	print("Euler: ",np.sum(XEerr2)**.5)
	print("RK2: ",np.sum(X2err2)**.5)
	print("Rk4: ",np.sum(X4err2)**.5)
	graph(XE,XEerr2,XE,X2err2,XE,X4err2,"Errors^2",'y','x(y) err^2',log=True)
	#graph(YE,XEerr2,Y2,X2err2,Y4,X4err2,"Errors^2",'y','x(y) err^2')


	return
if __name__ == "__main__":
	XE,YE,Y2,Y4 = Q2A()
	Q2B(XE,YE,Y2,Y4)
	Q2C(XE,YE,Y2,Y4)
