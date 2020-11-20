import numpy as np
from math import *
np.set_printoptions(suppress=True, precision=6)

#P36: Page 245 Problem 13
#Hint: See the predefined functions in 
#"module triangleQuad" provided in the book.
def triangleQuad(f,xc,yc):
  #area coordinates
	alpha = np.array([[1/3,1/3,1/3], \
										[.2, .2, .6], \
										[.6,.2,.2],\
										[.2,.6,.2]])
	W = np.array([-27,25,25,25])/48
	#Cubic weights
	x = np.dot(alpha,xc)#Convert to new coord basis
	y = np.dot(alpha,yc)
	#Find det for area coords?
	A = ( xc[1]*yc[2] - xc[2]*yc[1] \
			- xc[0]*yc[2] + xc[2]*yc[0] \
			+ xc[0]*yc[1] - xc[1]*yc[0])/2.0
	sum = 0.0
	for i in range(4):
		sum += W[i]*f(x[i],y[i])
	return A*sum

def P36():
	def f(x,y):
		return (1-x)*(y-x)*y
	xc = np.array([0,1,1.])
	yc = np.array([0,0,1.])
	print("+-----+")
	print("| P36 |")
	print("+-----+")
	print(triangleQuad(f,xc,yc))
	
def triTest():
	def f(x,y):
		return (x**2 + y**2)/2.0 \
					-(x**3 - 3.0*x*y**2)/6.0 \
					-2.0/3.0
	xCorner = np.array([ -1.0, 2.0, -1.0])
	yCorner = np.array([ -sqrt(3.0), 0.0, sqrt(3.0)])
	print("Integral =",triangleQuad(f,xCorner,yCorner))
	#-1.55884572681

#P37: Page 245 Problem 11
#Hint: See the predefined functions in "module gaussQuad2" 
#provided in the book.
from hw5 import gaussNodes

def gaussQuad2(f,x,y,m):
	def jac(x,y,s,t):#Find jacobian for change of basis
		J = np.zeros((2,2,))
		mt, pt, ms, ps = 1.0-t, 1.0+t, 1.0-s, 1.0+s
		J[0,0] = mt*(-x[0]+x[1]) + pt*(x[2]-x[3])
		J[0,1] = mt*(-y[0]+y[1]) + pt*(y[2]-y[3])
		J[1,0] = ms*(-x[0]+x[3]) + ps*(-x[1]+x[2])
		J[1,1] = ms*(-y[0]+y[3]) + ps*(-y[1]+y[2])
		#|J|
		return (J[0,0]*J[1,1]-J[0,1]*J[1,0])/16.0
	def map(x,y,s,t):
		N=np.zeros(4)
		#map coords to quadrilateral 
		mt, pt, ms, ps = 1.0-t, 1.0+t, 1.0-s, 1.0+s
		N[0] = ms*mt/4.0
		N[1] = ps*mt/4.0
		N[2] = ps*pt/4.0
		N[3] = ms*pt/4.0
		xCoord = np.dot(N,x)
		yCoord = np.dot(N,y)
		return xCoord,yCoord
	s,A = gaussNodes(m)#Find nodes
	sum = 0.0
	for i in range(m):
		for j in range(m):
			xCoord,yCoord = map(x,y,s[i],s[j])
			#Do summation with adjusted coords
			sum += A[i]*A[j]*jac(x,y,s[i],s[j])*f(xCoord,yCoord)
	return sum

def P37():
	def f(x,y):
		return x*y*(2-x**2)*(2-x*y)
	x = np.array([-3,1,3,-1.])
	y = np.array([-2,-2,2,2.])
	print("+-----+")
	print("| P37 |")
	print("+-----+")
	print(gaussQuad2(f,x,y,20))

def gQtest():
	def f(x,y):
		return (x**2 + y**2)/2.0 \
						- (x**3 - 3.0*x*y**2)/6.0 \
						- 2.0/3.0
	x = np.array([-1.0,2.0,2.0,-1.0])
	y = np.array([-sqrt(3.0),0.0,0.0,sqrt(3.0)])
	m = 3#eval(input("Integration order ==> "))
	print("Integral =", gaussQuad2(f,x,y,m))#-1.55884572681

#P38: Page 263 Problem 7
#Hint: See the predefined functions in "module run_kut4" 
#provided in the book. Also, look at Example 7.4 to have an idea.
def integrate(F,x,y,xStop,h):
	def run_kut4(F,x,y,h):
	  #RK4 Equations
		K0 = h*F(x,y)
		K1 = h*F(x + h/2.0, y + K0/2.0)
		K2 = h*F(x + h/2.0, y + K1/2.0)
		K3 = h*F(x + h, y + K2)
		return (K0 + 2.0*K1 + 2.0*K2 + K3)/6.0
	X, Y = np.array(x),np.array(y)
	while x < xStop: #Iterate through range (dont overstep)
		h = min(h,xStop-x)
		y = y + run_kut4(F,x,y,h)#find new y
		x = x + h#find new x
		X = np.append(X, x)
		Y = np.vstack((Y, y))
	return X,Y

import matplotlib.pyplot as plt

def P38():
	print("+-----+")
	print("| P38 |")
	print("+-----+")
	def F(x,y):#tau,theta
		F = np.zeros(2)
		F[0] = y[1]
		F[1] = -1*sin(y[0])#-0.1*y[1]-x
		return F
	x, xStop = 0.0, 3*pi#Range for X
	y = np.array([1.0,0.0])#Initial Conditions
	h = 0.001#Step size
	X,Y = integrate(F,x,y,xStop,h)
	for i in range(10,len(X)):
		if Y[i-1,1] >= 0 and Y[i,1]<=0:#Check to see if completed
		    #one full period
			T = X[i]
			print("Period:",round(T/pi,6),"pi*sqrt(L/g)")
			break
	plt.plot(X,Y[:,0],X,Y[:,1])#plot
	plt.grid(True)
	plt.xlabel('tau'); plt.ylabel('theta')
	plt.xlim(x,T)
	plt.savefig("hw6p38.png")#show()

#P39: Page 264 Problem 8
#Hint: See the predefined functions in "module run_kut4"
# provided in the book.
def P39():
	print("+-----+")
	print("| P39 |")
	print("+-----+")
	def F(x,y):#t,y
		F = np.zeros(2)
		g,CD,m=-9.80665,0.2028,80
		F[0] = y[1]
		F[1] = g+CD/m*y[1]*y[1]
		return F
	x, xStop = 0.0, 100# 0.5
	y = np.array([5000,0.0])#[1, 0])
	h = 0.1
	X,Y = integrate(F,x,y,xStop,h)
	for i in range(10,len(X)):
		if Y[i-1,0] >= 0 and Y[i,0]<=0:
			T = X[i]
			print("Fall time:",round(T,6),"s")
			break
	plt.clf()
	plt.plot(X,Y[:,0])
	plt.grid(True)
	plt.xlabel('t'); plt.ylabel('y')
	plt.xlim(x,T)
	plt.ylim(0,5000)
	plt.savefig("hw6p39.png")#show()

#P40: Page 268 Problem 22
#You do not need to submit the plot. 
#The program should contain a function which will return 
#two arrays containing the values of time and currents 
#respectively, used for the plot. Hint: See the predefined 
#functions in "module run_kut4" provided in the book.
def P40():
	print("+-----+")
	print("| P40 |")
	print("+-----+")
	def F(x,y):#y=[q1,i1,q2,i2]
		F = np.zeros(4)
		E,R,L,C=9,0.25,1.2*(10**-3), 5*(10**-3)
		F[0] = y[1]#q1'=i1
		F[1] = (E-R*y[1]-(y[0]-y[2])/C)/L#q1''=i1'
		F[2] = y[3]#q2'=i2
		F[3] = (-R*y[3]-(y[2]-y[0])/C-y[2]/C)/L#q2''=i2'
		return F
	x, xStop = 0.0, 0.05# 0.5
	y = np.array([0.0,0.0,0.0,0.0])#[1, 0])
	h = 0.0005
	X,Y = integrate(F,x,y,xStop,h)

	plt.clf()
	plt.plot(X,Y[:,1],X,Y[:,3])
	plt.grid(True)
	plt.xlabel('t'); plt.ylabel('i')
	plt.xlim(x,xStop)
	plt.legend(["i1","i2"])

	plt.savefig("hw6p40.png")#show()
	print(np.hstack((X.reshape(len(X),1),Y[:,(1,3)])))
	print("s, A, A")

if __name__ == "__main__":
	P36()
	P37()
	P38()
	P39()
	P40()
