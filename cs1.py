from math import *
from numpy import *
import time
#Simulation parameters: Range for x, t, and step size for each
Nx, Nt, Dx, Dt = 101, 9000, 0.03, 0.9
#Conductivity, specf heat, density - From measurements on 
#  material
KAPPA, SPH, RHO = 237, 900., 2700. 
#Solving the PDE using forward differention allows us to simplify
#   and reduce variables
MU = KAPPA / (SPH * RHO) * Dt / ( Dx**2 )
#Make matrix for Temperature - 
# use T[:, 1] for time now, 
#   T[:, 0] for time 1 dx ago (previous step)
#Tpl is for storing T every k iterations to aid in 
#   visualization and graphing
T, Tpl = zeros((Nx, 2), float), zeros((Nx, 31), float)
#LINE0: Insert intial and boundary conditions here
#Set the entire bar to 100c
T = T+100
#Set the ends of the bar to 0c
T[0,:], T[-1,:] = 0, 0

#Iterator for storing T in Tpl
m = 1
#Iterate through desired time range
t0 = time.time()
for t in range(1, Nt):
  #Iterate through the length of the bar but don't touch the ends
  # that are 0c
    for ix in range(1, Nx-1):
        #LINE1: Write update rule here
        #T(x,t+dx) = T(x, t) + MU * ( T(x+dx, t) + T(x-dx, t) 
        #   - 2*T(x,t) )
        T[ix, 1] = T[ix, 0] + MU  * ( T[ix+1, 0] + T[ix-1, 0] - \
           2*T[ix,0] )
    #Store T into Tpl so we can plot the time dependence
    if t%300 == 0 or t == 1:
        for ix in range(1, Nx - 1, 2): Tpl[ix, m] = T[ix, 1]
        print(m)
        m = m + 1
    #Update T so T(x,t)=T(x,t+dx) for next iteration
    for ix in range(1, Nx - 1): T[ix, 0] = T[ix, 1]
print("time:",time.time()-t0)

import matplotlib.pyplot as p
from mpl_toolkits.mplot3d import Axes3D

#Create axes for plots, use every other x
x, y = list(range(1, Nx - 1, 2)), list(range(1, 30))
#Create meshgrid for plotting
X, Y = meshgrid(x, y)

#Find corresponding T values in the mesh grid
def functz(Tpl):
    Z = Tpl[X, Y]
    return Z

Z = functz(Tpl)
#Create figure and plot X Y Z
fig = p.figure()
ax = Axes3D(fig)
ax.plot_wireframe(Y, X, Z, color = 'r')
ax.set_xlabel('time')
ax.set_ylabel('Position')
ax.set_zlabel('Temperature')
fig.suptitle("Heat Equation: Finite Difference")
p.savefig("figures/HE_finite_difference.png")
p.show()


#Create the range for x and t (needs to be symmetric?)
Max, n, m = 51, 50, 50
#Create vector for tridiagonal matrix a is diagonal below main in
#   A, b is standard B, 
#c is diag above main diag in A, d  is main diag in A
#T* is for storing and * is for operations
Ta, Tb, Tc, Td = zeros((Max),float), zeros((Max),float), \
  zeros((Max),float), zeros((Max),float)
a,   b,  c,  d = zeros((Max),float), zeros((Max),float), \
  zeros((Max),float), zeros((Max),float)
x, t = zeros((Max),float), zeros((Max,Max),float)

def Tridiag(a, b, c, d, Ta, Tb, Tc, Td, x, n):
    Max = 51
    #Create array for solved upper triangular values h and p
    #h will be primary diag of A, p is the new B
    h, p = zeros((Max),float), zeros((Max),float)
    for i in range(1, n+1):
      #reallocate data
        a[i], b[i], c[i], d[i] = Ta[i], Tb[i], Tc[i], Td[i]
    #solve for h1 and p1
    h[1], p[1] = c[1]/d[1], b[1]/d[1]
    for i in range(2, n+1):
        #LINE0: Put in equation for solving for h[i]
        #recursively find h_i from h1 to hn
        h[i] = c[i] / ( d[i] - a[i] * h[i-1] )
        #LINE1: Put in equation for solving for p[i]
        #recursively find p_i from p1 to pn
        p[i] = ( b[i] - a[i] * p[i-1] ) / ( d[i] - a[i] * h[i-1] )
    #we know last row in A has only one value of 1 in col N so 
    #   solve for Xn=pn
    x[n] = p[n]
    #LINE2: Put in equation for solving for x[i]
    #for i in range(1,n): x[i] = p[i] - h[i]*x[i-1]
    #use backwards substitution fron x[n-1] to x[1] knowing x[n] 
    #   is p[n]
    #each row has 1 on main diagonal and hi on main diag+1 => 
    #   x[i]+h[i]*x[i+1]=p[i]
    for i in range(n-1,0,-1): x[i] = p[i] - h[i]*x[i+1]

#create physical device parameters
width, height, ct= 1., .1, 1.
#for i in range(0, n): t[i,0] = 0.
#for i in range(1, m): t[0,i] = 0.
#h = dx, k = dt  
h, k = width / ( n - 1 ), height / ( m - 1 )
#r = eta = Conductivity*dt / ( SpecfHeat * density * dx^2 )
#r = KAPPA * k / (SPH * RHO * h**2 ) 
r = ct**2 * k / ( h**2 )

#Set the edges to 0c
for j in range(1,m+1):
    t[1,j], t[n,j] = 0., 0. #Boundary Condition
#LINE3:t[i][i] = Initial Condition = 100c
for i in range(2,   n): t[i,1]=100
#Set the main diagonal to 2/r + 2 => From Crank-nicolson PDE sln
for i in range(1, n+1): Td[i] = 2. + 2./r
#Make corners 1
Td[1], Td[n] = 1., 1.
#Set the off diagonals = -1
for i in range(1, n): Ta[i], Tc[i] = -1., -1.
#make the off diagonals that don't exist equal to 0, edges BC
Ta[n-1], Tc[1], Tb[1], Tb[n]  = 0., 0., 0., 0.
print("I'm working, wait for fig while I count to 50")

t1 = time.time()
#Iterate through matrix , dont touch BC
for j in range(2, m+1):
    print(j)
    for i in range(2, n):
      #reset Tb for next iteration
        Tb[i] = t[i-1,j-1] + t[i+1,j-1] + (2/r - 2) * t[i,j-1]
    #solve matrix
    Tridiag(a, b, c, d, Ta, Tb, Tc, Td, x, n)
    #store temp then into temp now
    for i in range(1, n+1): t[i,j] = x[i]
print("Finished")
print("time:",time.time()-t1)
#set up plotting matrix
x, y = list(range(1,m+1)), list(range(1, n+1))
X, Y = meshgrid(x, y)

Z = functz(t)
print(amin(Z))#Getting value below 0! Error with crank-nicolson 
#-> does average out though
fig = p.figure()
ax = Axes3D(fig)
ax.plot_wireframe(Y, X, Z, color = 'r')
ax.set_ylabel('Position')
ax.set_xlabel('time')
ax.set_zlabel('Temperature')
fig.suptitle("Heat Equation: Crank-Nicolson")
p.savefig("figures/HE_crank_nicolson.png")
p.show()

