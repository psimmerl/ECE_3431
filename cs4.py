from numpy import zeros, min, max, dot, \
    sum, linspace, arange, ones
from math import floor, ceil, exp, log, sqrt, pi, sin
import matplotlib.pylab as p
from mpl_toolkits.mplot3d import Axes3D

def plot(V, title, save_name, elev=20., azim=30.):
    #Get X,Y range
    NXmax, NYmax = V.shape
    x = range(0, NXmax-1, 1)
    y = range(0, NYmax-1, 1)
    X, Y = p.meshgrid(x,y)#Make grid
    def functz(V):#Make Z for plotting
        z = V[X,Y]
        return z
    Z = functz(V)
    fig = p.figure()#Make new fig
    ax = Axes3D(fig)
    #Plot fewer lines
    xstride, ystride = floor(max(x)/20), floor(max(y)/20)
    ax.plot_wireframe(X, Y, Z, rstride=xstride, \
        cstride=ystride, color = 'r')
    ax.set_xlabel('X'); ax.set_ylabel('Y')
    ax.set_zlabel('Potential, V(x,y)')
    ax.set_xlim(0,max(x)); ax.set_ylim(0,max(y))
    ax.set_title(title)
    ax.contourf(X,Y,Z,11,stride=25,offset=min(Z))
    ax.view_init(elev=elev, azim=azim)
    p.savefig('figures/%s.png'%save_name,bbox_inches='tight')
    print("Finished %s (%s.png)" % (title, save_name))

def poisson(V, rho=None, Niter=10000, tol=1.e-9):
    print("Initializing")
    NXmax, NYmax = V.shape#Get shape
    if rho is None:#make rho if nothing is passed
        rho = zeros(V.shape, float)
    print("Working hard, wait for figure while I count to %d" % \
        (floor((Niter-1)/10)*10))
    Vp = V.copy()
    for iter in range(Niter):#repeat many times
        if iter%10 == 0: print(iter)
        for i in range(1,NXmax-2):#x axis, skip BC
            for j in range(1,NYmax-2):#y axis, skip BC
                if rho[i,j] != 0:#If rho!=0 use rho eqn
                    V[i,j] = rho[i,j]
                elif abs(V[i,j]) != 100:#don't modify if IC
                    V[i,j] = (V[i+1,j]+V[i-1,j]+\
                        V[i,j+1]+V[i,j-1])/4#rho=0 eqn
        if sum((V-Vp)**2)/V.size < tol:#Check to see in tol
            print(iter)
            break
        Vp = V.copy()#make copy for tol checking

def gaussian(k, A=100, mu=0, sig=1):#Make gaussian change dist
    for i in range(len(k)):
        k[i] = A*exp(-((i-len(k)/2-mu+1/2)/sig)**2)
    return k

def parab(k, A=100):#Make parabolic charge dist
    for i in range(len(k)):
        k[i] = -i*(i-len(k)+1)/(len(k)/2-1/2)**2*A
    return k

def linear(k, A=100):#Make triangular charge dist
    for i in range(len(k)):
        k[i] = A*i/(len(k)/2)
        if i > len(k)/2:
            k[i] = A*(len(k)-i-1)/(len(k)/2-1/2)
        else:
            k[i] = A*i/(len(k)/2-1/2)
    return k

if __name__ == "__main__":
    print("Case Study 4")
    Nmax = 40#range of plot
    V = zeros((Nmax, Nmax), float)#make all 0, also BC=0
    V[0,1:(Nmax-1)] = 100 #Volts, IC-> could use rho but won't
    poisson(V)
    plot(V, "Single 100V wire on y axis", "cs4_1", azim=-45.)
    
    ang,el=40,20#set viewing angle
    L, w = 150, 60#For linear regime => d<<W->d~w/10
    d = 15
    xlb, xub, yp1, yp2 = floor((L-w)/2), ceil((L+w)/2), \
        floor((L-d)/2), ceil((L+d)/2)#set pos for plates
    V = zeros((L, L), float)
    V[xlb:xub, yp1] =  100 #Volts, IC-> could use rho but won't
    V[xlb:xub, yp2] = -100 #Volts, IC-> could use rho but won't
    poisson(V)
    plot(V, "Parallel Plate Capacitor", "cs4_2", elev=el, azim=ang)
   
    A, L = 100, 100 #set amplitude and range
    s = 25/sqrt(-log(0.01))#sigma for gauss, 0.01*MAX at ends
    V = zeros((L, L), float)#make initial grid+BC
    rho = zeros(V.shape,float)#Make rho
    rho[25:76,37] = gaussian(rho[25:76,37], A, 0, s)#set plate1
    rho[25:76,63] = -rho[25:76,37]#set plate2
    poisson(V,rho)
    plot(V, """Parallel Plate Capacitor w/ 
        Fitted Gaussian Charge Density""", "cs4_3", elev=el,azim=ang)

    V = zeros((L, L), float)#make initial grid+BC
    rho = zeros(V.shape,float)#Make rho
    rho[25:76,37] = gaussian(rho[25:76,37], A, 0, 2*s)#set plate1
    rho[25:76,63] = -rho[25:76,37]#set plate2
    poisson(V,rho)
    plot(V, """Parallel Plate Capacitor w/ 
        Large Gaussian Charge Density""", "cs4_4", elev=el,azim=ang)

    V = zeros((L, L), float)#make initial grid+BC
    rho = zeros(V.shape,float)#Make rho
    rho[25:76,37] = parab(rho[25:76,37], A)#set plate1
    rho[25:76,63] = -rho[25:76,37]#set plate2
    poisson(V,rho)
    plot(V, """Parallel Plate Capacitor w/ 
        Parabolic Charge Density""", "cs4_5", elev=el, azim=ang)

    V = zeros((L, L), float)#make initial grid+BC
    rho = zeros(V.shape,float)#Make rho
    rho[25:76,37] = linear(rho[25:76,37], A)#set plate1
    rho[25:76,63] = -rho[25:76,37]#set plate2
    poisson(V,rho)
    plot(V,"""Parallel Plate Capacitor w/ 
        Triangular Charge Density""", "cs4_6", elev=el,azim=ang)

    #Make plot of charge densities
    cs = ones(51,float)*100#Normal parallel plate
    fg = gaussian(zeros(51,float), A, 0, s)#fitted gaus
    lg = gaussian(zeros(51,float), A, 0, 2*s)#large gaus
    pb = parab(zeros(51,float), A)#parab
    tg = linear(zeros(51,float), A)#triang
    fig = p.figure()#make fig
    ax = p.gca()
    p.plot(arange(25,76),cs,arange(25,76),fg,arange(25,76),\
        lg,arange(25,76),pb,arange(25,76),tg)
    ax.set_xlabel('X'); ax.set_ylabel('Y')
    ax.set_xlim(25,75); ax.set_ylim(0,100.5)
    p.xticks(arange(25, 76, 5)); p.yticks(arange(0, 101, 10))
    ax.set_title("Charge Distributions")
    p.grid(True)
    p.legend(["Constant","Fitted Gaussian", "Large Gaussian", \
        "Parabolic", "Triangular"])
    p.savefig('figures/cs4_cd.png', bbox_inches='tight')

    #Tried sep of variables just to see soln
    L = 40
    sC = pi/L
    C1 = 100/(1-exp(-2*sC*L))
    V = zeros((L, L), float)
    for x in range(L):
        for y in range(L):
            V[x,y]=(C1*exp(-sC*x)+(100-C1)*exp(sC*x))*sin(sC*y)
    plot(V, """Analytic Sep
     Vars for Single 100V 
     wire on y axis""", "cs4_1a", azim=-45.)

    


    
