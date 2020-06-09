import scipy as sp
import numpy as np
import pylab as pl
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as spi

G= 6.67*10**(-11) #define G - constant of gravity
M = 6.4*10**(23)  #define M - mass of Mars
m = 260     #define m - mass of satellite
R = 3.4*10**(6)   #define mars’ radius
mv = 24100  #mars’ velocity
marsmoving = input(str("Do you want to simulate mars moving? Enter y for yes, enter n for no:")) #allows user to
def f(x,t): #called with the variables to be integrated
    if marsmoving == "y":
        xx=x[0] - mv*t #updating satellites’ position relative to mars
    elif marsmoving == "n":
        xx=x[0] # updates position of satellite relative to origin
    vx=x[1] # x component of velocity
    yy=x[2] # y component of distance
    vy=x[3] # y component of velocity
    r = (xx**2 +yy**2)**0.5   #radial distance from centre of Mars
    if r < R: #dont plot - zero values for [xx,vx,yy,vy]
        if marsmoving == "y": #keeps satellite moving with mars’ velocity when it crashes
            return[mv,0,0,0]
        elif marsmoving == "n": #keeps satellite stationary when crashes on mars
            return[0,0,0,0]
    ax=-xx*(G*M)/(r**3)  #acceleration in x direction
    ay=-yy*(G*M)/(r**3) #acceleration in y direction
    return [vx,ax,vy,ay] #actually returns [xx,vx,yy,vy] - integrates returned values
timelimit = int(input("Please enter the time limit you would like to run the simulation:")) #user option for time sca
initialv = int(input("Enter initial velocity:")) #user option for initial velocity
finalv = int(input("Enter final velocity:")) #user option for final velocity
t=sp.linspace(0,timelimit,4000) # solving every t per second
Ang = [] # empty list for angle deviation
vlist = [] # empty list for velocites to be used for angle deviation
dist =[]  # for minimum dist

for i in np.arange(initialv, finalv,10):
    if marsmoving == "y":
        IC=[0,mv+i,4*R,i] # set initial conditions [x,vx,y,vy] for satellite when mars is moving
    elif marsmoving =="n":
        IC=[7*R,-i,3*R,i] # set initial conditions [x,vx,y,vy] for satellite when mars is stationary
    soln=spi.odeint(f,IC,t) # integrates vx,vy,ax,ay
    x=soln[:,0] #takes xx solution and places it in first column of array
    vx=soln[:,1] #takes vx solution and places it in second column of array
    y=soln[:,2] #takes yy solution and places it in third column of array
    vy=soln[:,3] #takes vy solution and places it in fourth column of array
    c = pl.get_cmap('Reds')((i - initialv)/finalv) # colouring different velocities, second brackets give rgb co
    dotproduct = (vx[0]*vx[-1]+vy[0]*vy[-1]) / (((vx[0]**2+vy[0]**2)**0.5)*((vx[-1]**2+vy[-1]**2)**0.5)) #dot pr
    O = 180*np.arccos(dotproduct)/np.pi # angle deviation for specific initial conditions
    Ang.append(O) #array of angle deviations
    v = (vx[0]**2+vy[0]**2)**0.5   # initial velocities
    vlist.append(v) #array of initial velocities to be used for angle deviation plot
    rdist =(x**2+y**2)**0.5 # for use in closest aproach
    if rdist[i] > R:
        dist.append(rdist)
    if marsmoving == 'n':
        plt.plot(x,y,color=c) # Plot the x and y components of satellite for multpile trajctories in the initial
        if finalv <650: #this was the finalv that gave close aproach
            zz = ", Closest distance:" + str(int(np.amin(dist))-int(R)) #calculates the minimum distance from th
        else:
            zz = ""
yM = [] #empty array for y position of mars
xM = [] #empty array for x position of mars
for i in range(len(t)): #for loop that plots mars’s position at same time as the satelittes motion
    yM.append(0) #keep adding zero’s to mars’ y position to keep mars on the x-axis so moves straight line
    xM.append(mv*t[i]) # appends the updated version of mars’ position as time moves forward
plt.figure(1)
if marsmoving == "y":
    plt.plot(x,y)
    plt.plot(xM,yM, 'ro') # plots mars moving
elif marsmoving == "n":
    mars=plt.Circle((0,0),R,color='r',label="Mars") # creates circle radius of mars
    plt.gcf().gca().add_artist(mars) # plots the cirlce
plt.xlabel("X position of Satellite",size=15)
plt.ylabel("Y position of Satellite",size=15)
if marsmoving == "n":
    if finalv < 600:                        #
        plt.xlim(-1e7,3e7)                  #
        plt.ylim(-1e7,3e7)                  #
    elif finalv >= 600 and finalv < 1200: ## setting the time scale for plots based on final velocity as this
        plt.xlim(-5e7,8e7)                  #
        plt.ylim(-5e7,8e7)                  #
    elif finalv > 1200:                     #
        plt.xlim(-2e8,2e8)                  #
        plt.ylim(-2e8,2e8)                  #
if marsmoving == "n":
    plt.title('Satellite moving around Mars' + zz, size=17)
else:
    plt.title('Satellite moving around Mars', size=17)
fig = plt.gcf() #takes figure plotted and lets edit
fig.set_size_inches((8, 8))
fig.savefig("Satellite’s trajcetory")
plt.show()
plt.figure(2) ##################Angle Deviation------------------------------------------
plt.plot(vlist,Ang, color = 'c')
plt.xlabel("Initial launch velocity (m/s)",size=15)
plt.ylabel("Angle deviation (Degrees)",size=15)
plt.title("Angle deviation for satellite",size=17)
fig = pl.gcf() #takes figure plotted and lets edit
plt.savefig("Angle deviation")
fig.set_size_inches((8, 8)) #sets size of figure
plt.figure(3)##################Energy Plots ---------------------------------------------
EK = 0.5*m*(vx**2+vy**2) # kinetic energy formula
r = (x**2+y**2)**0.5 #radial distance
PE = -G*M*m/r  # potential energy
plt.plot(t,EK, label = "EK") #plot kinetic
plt.plot(t,PE, label = "PE") #plot potential
plt.plot(t,EK+PE, label ="EK + PE") # plot total enegry
plt.title("Energy vs. Time for captured satellite ahead of Mars",size=13)
matplotlib.rcParams.update({'font.size': 13})
plt.xlabel("Time (Seconds)",size=15)
plt.ylabel("Energy (Joules)",size=15)
plt.legend()
fig = plt.gcf() #takes figure plotted and lets edit
fig.set_size_inches((8, 8))
plt.savefig("Energy Plot")
plt.show()
