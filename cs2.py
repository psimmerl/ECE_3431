#Case Study 2
from math import *
import numpy as np

def LUdecomp3(c,d,e):
  n = len(d)
  for k in range(1,n):
    lam = c[k-1]/d[k-1]
    d[k] -= lam*e[k-1]
    c[k-1] = lam
  return c,d,e
def LUsolve3(c,d,e,b):
  n = len(d)
  for k in range(1,n):
    b[k] -= c[k-1]*b[k-1]
  b[n-1] = b[n-1]/d[n-1]
  for k in range(n-2, -1, -1):
    b[k] = (b[k]-e[k]*b[k+1])/d[k]
  return b

def curvatures(xData, yData):
  n = len(xData) - 1
  #Create Tridiag matrix with d on main diag
  #c on upper diag, e on lower diag and k as B
  #see eq 12
  c, d = np.zeros(n), np.ones(n+1)
  e, k = np.zeros(n), np.zeros(n+1)
  c[0:n-1] = xData[0:n-1] - xData[1:n]
  d[1:n] = 2.0*(xData[0:n-1] - xData[2:n+1])
  e[1:n] = xData[1:n] - xData[2:n+1]
  k[1:n] = 6.0*(yData[0:n-1] - yData[1:n]) \
            /(xData[0:n-1] - xData[1:n])  \
          -6.0*(yData[1:n] - yData[2:n+1]) \
            /(xData[1:n] - xData[2:n+1]) 
  #Decompose tridiag and solve
  LUdecomp3(c,d,e)
  LUsolve3(c,d,e,k)
  return k

def evalSpline(xData,yData,k,x):
  #Determine which f_{i,i+1}(x) is used
  def findSegment(xData,x):
    iLeft = 0
    iRight = len(xData)-1
    while True:#Binary bisection
      #If stuck between 2 f, use left f
      if (iRight-iLeft) <= 1: return iLeft
      #choose middle data point
      i = int((iRight-iLeft)/2+iLeft)
      #check to see if x is less than data point
      if x < xData[i]: iRight = i
      else: iLeft = i

  #Find which f
  i = findSegment(xData,x)
  #Make math easier to read
  h = xData[i] - xData[i+1]
  #Use equation 11 to evaluate spline
  y = ((x-xData[i+1])**3/h - (x - xData[i+1])*h)*k[i]/6.0 \
      -((x-xData[i])**3/h - (x - xData[i])*h)*k[i+1]/6.0   \
      + (yData[i]*(x-xData[i+1]) \
      - yData[i+1]*(x-xData[i]))/h
  return y

import matplotlib.pyplot as plt
#Plotting spline and data helper function
def myPlot(X, Y, xs, ys, title, fname, xlabel='X', ylabel='Y'):
  plt.clf()
  plt.plot(xs,ys)#Plot spline
  plt.plot(X,Y,'rx')#Plot discrete data
  plt.grid()
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  plt.savefig(fname)
  #plt.show()

print("+----------------+")
print("|   Question 1   +")
print("+----------------+")
X = np.array([1,2,3,4,5])
Y = np.array([0,1,0,1,0])
k = curvatures(X,Y)#Find spline
print(k)
print(evalSpline(X,Y,k,1.5))
#Make test x from 0.1 to 5 on .1 intervals
xs = np.arange(0.1, 5.1, 0.1)
ys = np.array([0]*len(xs),dtype=float)
for i in range(len(xs)):#Eval each x with spline
  ys[i] = evalSpline(X,Y,k,xs[i])
asort = ys.argsort()#Sort x and y based on y
xas = xs[asort]
yas = ys[asort]
for i in range(len(xas)):#iterate through x
  tol = 1e-6
  #if y is within tol of surrounding data, print x,y
  if i > 1 and abs(yas[i-1]-yas[i])<tol: print(xas[i], yas[i])
  elif i<(len(xas)-1) and abs(yas[i+1]-yas[i])<tol: print("\n",xas[i], yas[i])
#Make plot
myPlot(X,Y,xs,ys,'Case Study 2 Problem 1','cs2p1.png')

print("+----------------+")
print("|   Question 2   +")
print("+----------------+")
X = np.array([0,1,2.])
Y = np.array([0,2,1])
k = curvatures(X,Y)#find spline
print(k)
tol,di,found = 1e-6,0.01, False
while di > tol:#if not finding data, try smaller points
  for i in np.arange(0, 1+di, di):#keep 0<=i<=1
    j = evalSpline(X,Y,k,i)#eval i for j
    if j <= 2 and j >= 1:#if j in range 1<=j<=2
      i2 = evalSpline(X,Y,k,j)#eval j for i2
      if abs(i2-i)<tol: #if eval j is within tol of i
        print(i,j)
        found = True#don't need more data
        break
  if found or di < tol: break#Finish if found a point
                              # or too small
  else: di = di/10#If didn't find, try more data

#Make spline for plotting
xs = np.linspace(min(X),max(X),100)
ys = np.array([0]*len(xs),dtype=float)
for i in range(len(xs)):
  ys[i] = evalSpline(X,Y,k,xs[i])
myPlot(X,Y,xs,ys,'Case Study 2 Problem 2','cs2p2.png')

print("+----------------+")
print("|   Question 3   +")
print("+----------------+")
#Make function for easy eval
def evalQ3(x):
  X = np.array([ 1, 2, 3,4, 5.])
  Y = np.array([13,15,12,9,13.])
  k = curvatures(X,Y)
  return evalSpline(X,Y,k,x)

print(evalQ3(3.4))
xs = np.linspace(min(X),max(X),100)
ys = np.array([0]*len(xs),dtype=float)
for i in range(len(xs)):
  ys[i] = evalQ3(xs[i])
myPlot(X,Y,xs,ys,'Case Study 2 Problem 3','cs2p3.png')
#plt.show()



print("+----------------+")
print("|   Question 4   +")
print("+----------------+")
#FIND Y(X)
X = np.array([  0.2,  0.4,  0.6,   0.8,   1.0])
Y = np.array([1.150,0.855,0.377,-0.266,-1.049])
#Resort data so we can find spline Y(x)
asort = Y.argsort()
X = X[asort]
Y = Y[asort]
k = curvatures(Y,X)
print(k)

#Find zeros (roots) of Y(x)
from hw3 import multisearch_ridder
def f(y): return evalSpline(Y,X,k,y)
multisearch_ridder(f,1,3)

#Plot spline with roots
ys = np.linspace(min(Y),max(Y)*1.9,100)
xs = np.array([0]*len(ys),dtype=float)
for i in range(len(ys)):
  xs[i] = evalSpline(Y,X,k,ys[i])
myPlot(Y,X,ys,xs,'Case Study 2 Problem 4','cs2p4.png','Y','X')

#plot spline show behavoior is unpredicable outide of data range
ys2 = np.linspace(floor(min(Y)*3),ceil(max(Y)*3),200)
xs2 = np.array([0]*len(ys2),dtype=float)
for i in range(len(ys2)):
  xs2[i] = evalSpline(Y,X,k,ys2[i])
myPlot(Y,X,ys2,xs2,'Case Study 2 Problem 4 Wide View','cs2p4_wide.png','Y','X')
plt.xlim(min(ys2),max(ys2))
plt.savefig('cs2p4_wide.png')



