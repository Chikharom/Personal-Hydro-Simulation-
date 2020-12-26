import numpy as np
import matplotlib.pyplot as plt
import math
n = 81
def generate_grid(size):
    '''Function to generate the grid on which we simulate the fluid. it also specifies
    which points are considered to be an obstacle (here a diamond)'''
    x = np.linspace(-4,4,size)
    y = np.linspace(-4,4,size)
    last_point = x[-1]
    xx, yy = np.meshgrid(x,y)
    fixed = np.zeros(shape=xx.shape)# set to 1 if the point is fixed
    #The following loop is specific to a diamond shaped obstacle in the middle
    #Use of 4 lines to define the diamond (can easily be modified for all other shapes)
    for i in range(len(xx)):
        for j in range(len(xx[0])):

            x_current = xx[i][j]
            y_current = yy[i][j]

            if xx[i][j] == -4 or x_current == last_point or y_current == -4 or y_current == last_point:

                fixed[i][j] = 1
            elif x_current >= -2 and x_current <= 0 and y_current <= x_current +2 and  y_current >= -x_current -2:

                fixed[i][j] = 1
            elif x_current >= 0 and x_current <= 2 and y_current <= -x_current +2 and  y_current >= x_current -2:

                fixed[i][j] = 1

    return xx,yy,fixed

xx,yy,fixed = generate_grid(n)

def initialize_grid(xx,yy):
    ''' Generated a potential function through the use of laplace's equation in
    a stationary flow. We leave each point at value 0 is it's part of the obstacle'''
    flow = np.zeros(shape=xx.shape)
    for i in range(len(xx)):
        for j in range(len(xx[0])):
            x_current = xx[i][j]
            y_current = yy[i][j]
            #checks obstacle
            if x_current >= -2 and x_current <= 0 and y_current <= x_current +2 and  y_current >= -x_current -2:
                pass
            elif x_current >= 0 and x_current <= 2 and y_current <= -x_current +2 and  y_current >= x_current -2:
                pass
            else :
                #Boundary conditions can be varied but are fixed here
                flow[i][j] = yy[i][j]
    return flow
alpha = 0.001
flow = initialize_grid(xx,yy)
def update(flow,fixed):
    ''' Function that updates the value of the potential at each point based on the
    value of the function evaluated at its neighbors. (Solving laplace's equation)'''
    new = np.copy(flow)#Ancien flow
    for i in range(len(flow)):
        for j in range(len(flow[0])):
            if fixed[i][j] == 1:#Fixed point

                pass
            else :
                #Give it the average of its neighbors
                new[i][j]=(flow[i+1][j]+flow[i-1][j]+flow[i][j+1]+flow[i][j-1])/4
    return new
first = update(flow,fixed)-flow
max = first.max()
min = first.min()
def check_max(a,b):
    if a > b:
        return a
    else :
        return b
current_max = check_max(max,abs(min))
def get_speed(potential):
    ''' Function to obtain speed from potentiel considering that the curl of the
    potential will give us the vector field here '''
    N = len(potential)
    speeds = np.zeros(potential.shape)#Magnitude of velocities
    speeds_x = np.zeros(potential.shape)#x velocities
    speeds_y = np.zeros(potential.shape)#y velocities
    for i in range(1,len(potential)):
        for j in range(1,len(potential[0])):
            #Finite difference to evaluate the necessary derivatives
            speeds_x[i][j] = (potential[i][j]-potential[i-1][j]) / (8/N)
            speeds_y[i][j] = -(potential[i][j]-potential[i][j-1]) / (8/N)
            speeds[i][j] = speeds_x[i][j]**2 + speeds_y[i][j]**2

    return speeds, speeds_x, speeds_y

def get_pression(speeds):
    ''' Obtaining pressure from speed using the bernouli relation'''

    pression = np.zeros(speeds.shape)
    for i in range(len(pression)):
        for j in range(len(pression[0])):
            #Theoretical result
            pression[i][j] = (-1/2)*speeds[i][j]

    return pression
def draw(xx,yy,potential,pression):
    '''Simple function to plot speeds and pressure to get an idea of the
    behaviour of the fluid'''
    y1 = [0, 2]
    y2 = [2, 0]
    y3 = [-2, 0]
    y4 = [0, -2]
    x1=[-2, 0]
    x2=[0, 2]
    x3=[0, 2]
    x4=[-2, 0]
    plt.plot(x1, y1, 'mediumseagreen', linewidth = 1.9)
    plt.plot(x2, y2, 'mediumseagreen', linewidth = 1.9)
    plt.plot(x3, y3, 'mediumseagreen', linewidth = 1.9)
    plt.plot(x4, y4, 'mediumseagreen', linewidth = 1.9)
    plt.pcolor(xx, yy, potential)
    plt.contour(xx, yy, potential, levels = 25, colors =['mediumseagreen', 'mediumseagreen'], alpha = 1)
    plt.contourf(xx,yy, pression, levels = 50, colors = ['palegreen', 'darkolivegreen', 'paleturquoise'])
    plt.fill_between(x1, y1, color='grey')
    plt.fill_between(x2, y2, color='grey')
    plt.fill_between(x3, y3, color='grey')
    plt.fill_between(x4, y4, color='grey')
    plt.show()
#Lines used to the integral of the drag force
def first_line(x):
    return x + 2
def second_line(x):
    return -x-2
def third_line(x):
    return 2-x
def fourth_line(x):
    return x-2
def get_drag_force(xx,yy,flow,it,show=True):
    ''' Using normal vectors to each surface of our obstacle, this function calculates the drag
    force experienced by the obstacle using the pressure at each point. Here it is implemented for the
    diamond obstacle.'''
    first_vector = np.array([-1,1])/ math.sqrt(2)
    second_vector = np.array([-1,-1])/ math.sqrt(2)
    third_vector = np.array([1.1])/math.sqrt(2)
    fourth_vector = np.array([1,-1]) / math.sqrt(2)
    speeds,speeds_x,speeds_y = get_speed(flow)
    pression = get_pression(speeds)
    integrale = 0
    beta = (8 / n)*(1/2) #Distance max to the border of the obstacle
    dS = beta/2
    counter = 0
    for i in range(len(flow)):
        for j in range(len(flow[0])):
            #Dot product with the right normal vector if the point is close enough
            if xx[i][j] >= -2 and xx[i][j] <= 0:

                if abs(first_line(xx[i][j])-yy[i][j]) <= beta:

                    integrale += pression[i][j]*first_vector[0]
                    counter +=1
                if abs(second_line(xx[i][j])-yy[i][j]) <= beta:

                    integrale += pression[i][j]*second_vector[0]
                    counter +=1
            elif xx[i][j] > 0 and xx[i][j] <= 2 :

                if abs(third_line(xx[i][j])-yy[i][j]) <= beta:

                    integrale += pression[i][j]*third_vector[0]
                    counter +=1
                if abs(fourth_line(xx[i][j])-yy[i][j]) <= beta:

                    integrale += pression[i][j]*fourth_vector[0]
                    counter +=1
    #Multiply by the step taking dS
    integrale = integrale*dS

    return integrale

def check_max(a,b):
    if a >b :
        return a
    else :
        return b
def Iteration_number(alpha,flow):
    '''function to determine the number of iterations required to have a certain
    precision alpha'''
    i = 0#_Iterations
    past = np.copy(flow)# Premier potentiel
    flow = update(flow,fixed)#2ieme potentiel
    temp = flow-past#difference
    #max de difference
    current_max =check_max(temp.max(),abs(temp.min()))
    while(current_max > alpha):
        past = np.copy(flow)
        flow = update(flow,fixed)
        temp = flow-past
        current_max =check_max(temp.max(),abs(temp.min()))
        i +=1
    return i
def Iteration(Iter,flow):
    for it in range(1,Iter):
        if (it%100==0):
            get_drag_force(xx,yy,flow,it)
        flow = update(flow,fixed)
    return flow
Iter = Iteration_number(0.0001,flow)
Iter = 2000
flow = Iteration(Iter,flow)
speeds, speeds_x, speeds_y = get_speed(flow)
pression = get_pression(speeds)
draw(xx,yy,flow,pression)
#plt.savefig("Graph_with_"+str(Iter+"_Iterations")
#plt.contour(xx,yy,get_speed(flow),levels=20)
#plt.show()
