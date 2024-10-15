import sys
import rebound
import reboundx
import numpy as np
import pandas as pd
import math
import time
import matplotlib.pyplot as plt
from tqdm import tqdm
#### Some Numbers
R_sun = 2.254e-6 # in 0.01 parsecs = 1 a
TDR = 2.254e-4 # in 0.01 parsecs = 1 a
e_sun = 0.9997746
M_sun = 1e-6
M_BH = 1 # 10^6 M_sun BH
delta_epsilon = 44.36 # M_bh*(0.01 parsecs)^2/(Some Fucked up unit) or 1.9e13 Joules
N = 2 # Number of splits


###########################################
# EARLY ACCESS TOOLS
###########################################

def trial_function(x,k,a,mu):
    return a*(x-mu)**-k

def uniform(x,a):
    return np.ones(len(x))*a

def GenerateParameter(N,f,a_in,a_out,**kwargs):
    f_max = -math.inf
    f_min = math.inf
    steps= np.linspace(a_in,a_out,10000)
    sub_div = np.diff(steps)[0]
    func_points = f(steps,**kwargs)
    area = np.sum(func_points[:-1]+func_points[1:])/2*sub_div

    for i in steps:
        if(f(i,**kwargs) > f_max):
            f_max = f(i,**kwargs)
        if(f(i,**kwargs) < f_min):
            f_min = f(i,**kwargs)

    final = np.array([])
    for i in range(N):
        found = False
        while(found == False):
            a = np.random.uniform()*(a_out-a_in)+a_in
            y = (np.random.uniform()*f_max)/area
            if(y <= f(a,**kwargs)/area):
                final = np.append(final,a)
                found = True

    return final

def GenerateDistribution(N,f,x_min=0,x_max=1, **kwargs):
    edges = np.linspace(x_min,x_max,N+1)
    x = np.array([])
    for i in range(N):
        points = np.array([edges[i],edges[i+1]])
        x = np.append(x,np.mean(points))
    y = f(x, **kwargs)
    return x,y

##################################################
# COLLISION FUNCTIONS
##################################################

def vampire(sim_pointer, collided_particles_index):
    global N
    colors = [(1.,0.,0.),(0.,0.75,0.75),(0.75,0.,0.75),(0.75, 0.75, 0,),(0., 0., 0.),(0., 0., 1.),(0., 0.5, 0.)]

    print('collision')
    sim = sim_pointer.contents # retreive the standard simulation object
    ps = sim.particles # easy access to list of particles
    print(len(ps),N)

    if(ps[collided_particles_index.p1] == ps[0]):
        i = collided_particles_index.p1   # Note that p1 < p2 is not guaranteed.    
        j = collided_particles_index.p2 
    else:
        j = collided_particles_index.p1   
        i = collided_particles_index.p2 

    distance = np.sqrt((ps[i].x-ps[j].x)**2+(ps[i].y-ps[j].y)**2)
    # sim.N_active = len(ps)

    # This part is exciting! We can execute additional code during collisions now!
    op = rebound.OrbitPlot(sim, particles = [j], primary = 0, color=True)
    op1 = rebound.OrbitPlot(sim, fig=op.fig, ax = op.ax)
    op.ax.set_title("Particle {} and {} collision".format(j, i))
    op.ax.text(ps[i].x, ps[i].y, "1") 
    op.ax.text(ps[j].x, ps[j].y, "2");
    plt.show()
    op2 = rebound.OrbitPlot(sim, particles = [j], primary = 0, color=True)
    op2.ax.set_xlim(-1.5*distance,1.5*distance)
    op2.ax.set_ylim(-1.5*distance,1.5*distance)
    op2.ax.set_title("Particle {} and {} collision zoomed".format(j,i))
    op2.ax.text(ps[i].x, ps[i].y, "1") 
    op2.ax.text(ps[j].x, ps[j].y, "2");
    plt.show()
    # So we plot the scenario exactly at the timestep that the collision function is triggered
    
    ##################### Splitting ###########################
    num_splits = N

    epsilon_range = np
    
    reduced_mass = ps[j].m/num_splits
    x = ps[j].x
    y = ps[j].y
    z = ps[j].z
    vx = ps[j].vx
    vy = ps[j].vy
    vz = ps[j].vz
    v_vec = np.array([vx,vy,vz])
    v_mag = np.sqrt(vx**2+vy**2+vz**2)
    v_hat = v_vec/np.sqrt(vx**2+vy**2+vz**2)

    # E_dist,trash = GenerateDistribution(num_splits,uniform,x_min=-delta_epsilon,x_max=delta_epsilon,a=1)
    temp_delta_epsilon = (1/2)*(0.1*v_mag)**2*reduced_mass
    E_dist,trash = GenerateDistribution(num_splits,uniform,x_min=-temp_delta_epsilon,
                                        x_max=temp_delta_epsilon,a=1)

    E_sign = E_dist/np.abs(E_dist)
    E_dist = E_dist/E_sign
    v_dist = np.sqrt(2*E_dist/reduced_mass)
    v_dist = v_dist*E_sign
    # print(v_vec)
    # print(v_dist)
    v_new = np.zeros((num_splits,3))
    for k in range(num_splits):
        v_new[k] = v_hat*v_dist[k]+v_vec
    # print(v_new)
    # a = ps[j].a
    # e = ps[j].e
    # inc = ps[j].inc
    # omega = ps[j].omega
    # Omega = ps[j].Omega
    # M = ps[j].M

    plot_list = []
    for k in range(num_splits):
        # sim.add(m=reduced_mass,\
        #     a=(a*0.5*(i+1)),e=e,M=M,inc=inc,Omega=Omega,omega=omega,jacobi_masses=False)
        sim.add(m=reduced_mass,x=x,y=y,z=z,vx=v_new[k,0],
                vy=v_new[k,1],vz=v_new[k,2],jacobi_masses=False)
        plot_list.append(-1*(k+1))

    plot_list_array = np.array(plot_list)
    plot_list_array = np.sort(plot_list_array)
    
    orbits = sim.orbits()
    print(len(orbits))
    e = []
    for k,orb in enumerate(orbits):
        if((k == j) or (k >= len(ps)-len(plot_list))):
            print(f"particle {k}, e = {orb.e}",f"e_x = {orb.evec.x}",f"e_y = {orb.evec.y}")
        e.append(orb.e)
    
    e_sub = [e[x] for x in plot_list_array]
    e_sub = np.array(e_sub)
    unbound_mask = e_sub > 1
    unbound = plot_list_array[unbound_mask]
    bound = plot_list_array[~unbound_mask]
    unbound = [x.item() for x in plot_list_array[unbound_mask]]
    bound = [x.item() for x in plot_list_array[~unbound_mask]]
    print(bound, unbound)
    unbound_colors = []
    bound_colors = []
    for k in range(len(unbound)):
        unbound_colors.append(colors[0])
    for k in range(len(bound)):
        bound_colors.append(colors[1])

    #################Plotting######################################
    # op = rebound.OrbitPlot(sim, particles = [j], primary = 0, lw=2)
    # op1 = rebound.OrbitPlot(sim, particles = plot_list, primary = 0,\
    #                         color = True, fig=op.fig, ax = op.ax)
    # op.ax.set_title("Particle {} split".format(j))
    # op.ax.text(ps[j].x, ps[j].y, j) 
    # plt.show()
    
    op2 = rebound.OrbitPlot(sim, particles = unbound, primary = 0,\
                            color = unbound_colors)
    op3 = rebound.OrbitPlot(sim, particles = bound, primary = 0,\
                            color = bound_colors, fig = op2.fig, ax = op2.ax)
    op4 = rebound.OrbitPlot(sim, particles = [j], primary = 0, lw=2,\
                            fig=op2.fig, ax = op2.ax)
    op2.ax.set_title("Particle {} split zoomed".format(j))
    op2.ax.set_xlim(-1.5*distance,1.5*distance)
    op2.ax.set_ylim(-1.5*distance,1.5*distance)
    op2.ax.scatter([],[],color=colors[0],label="Unbound",marker='_')
    op2.ax.scatter([],[],color=colors[1],label="Bound",marker='_')
    op2.ax.text(ps[j].x, ps[j].y, j)
    op2.ax.legend()
    plt.show()

    op5 = rebound.OrbitPlot(sim, particles = bound, primary = 0,\
                            color = bound_colors)
    op6 = rebound.OrbitPlot(sim, particles = [j], primary = 0,\
                            fig = op5.fig, ax = op5.ax)
    op7 = rebound.OrbitPlot(sim, particles = unbound, primary = 0, lw=2,\
                            color = unbound_colors, fig=op5.fig, ax = op5.ax)
    op5.ax.set_title("Particle {} split outcome".format(j))
    op5.ax.scatter([],[],color=colors[0],label="Unbound",marker='_')
    op5.ax.scatter([],[],color=colors[1],label="Bound",marker='_')
    op5.ax.text(ps[j].x, ps[j].y, j)
    op5.ax.legend()
    plt.show()

    op8 = rebound.OrbitPlot(sim, particles = unbound, primary = 0,\
                            color = unbound_colors,lw = 10)
    op9 = rebound.OrbitPlot(sim, particles = bound, primary = 0,\
                            color = bound_colors, lw = 10, fig = op8.fig, ax = op8.ax)
    op10 = rebound.OrbitPlot(sim,primary = 0, fig = op8.fig, ax = op8.ax, alpha = 0.1)

    op8.ax.set_title("Particle {} split Unbound".format(j))
    op8.ax.scatter([],[],color=colors[0],label="Unbound",marker='_')
    op8.ax.scatter([],[],color=colors[1],label="Bound",marker='_')
    op8.ax.text(ps[j].x, ps[j].y, j)
    op8.ax.legend()
    plt.show()

    op11 = rebound.OrbitPlot(sim, primary = 0,alpha = 0.1)
    op12 = rebound.OrbitPlot(sim, particles = bound, primary = 0,\
                            color = bound_colors, lw = 10, fig = op11.fig, ax = op11.ax)
    op113 = rebound.OrbitPlot(sim,particles = unbound, primary = 0,\
                              color = unbound_colors, lw = 10, fig = op11.fig, ax = op11.ax)

    op11.ax.set_title("Particle {} split Totality".format(j))
    op11.ax.scatter([],[],color=colors[0],label="Unbound",marker='_')
    op11.ax.scatter([],[],color=colors[1],label="Bound",marker='_')
    op11.ax.text(ps[j].x, ps[j].y, j)
    op11.ax.legend()
    plt.show()


    sys.exit()
    # for i in range(num_splits):
    #     ps = sim.particles
    #     sim.remove(len(ps)-1)

    # sim.remove(plot_list)

    # sys.exit()
    if(i == 0):
        print("Return",i,j)
        return 2 # remove particle with index j
    elif(j == 0):
        print("return",i,j)
        return 1
    else:
        print("Return 0")
        return 0
    # return 0 # Don't remove either particle

###############################################################
# STARTING FUNCTION
###############################################################

def Initialize(massed=True,collisions=False,N_splits=2): 
    # massed is an argument for if you are loading a previously made sim or if it is one created from scratch
    # if loading sime, massed=False, if from scratch, massed=True. Also not necessary because it will be overwrittin if loading
    global N
    N = N_splits
    sim = rebound.Simulation()
    if(collisions):
        sim.collision = 'line'
        sim.collision_resolve = vampire
    if(massed):
        # sim.add(m=1,r=2/(1.37*10**4)**2)
        sim.add(m=1,r=TDR)
    return sim

##################################################################
# GENERATING A POPULATION OF PARTICLES IN DIFFERENT CONFIGURATIONS
##################################################################

def GenerateDisk(sim,N,a_range,massless=False):
    for i in tqdm(range(N)):
        if(massless==False):
            m = M_sun
            # r = 2*m/(1.37*10**4)**2 # c = 1.37*10^4 (radius units)/(time units) assuming
            # a period around central object is (2 * pi) = (1 year), (1 radius) = (1 au)
        else:
            m=0
            # r = 2*(1/332900)/(1.37*10**4)**2 # same as before but assuming mass of earth
        a = GenerateParameter(1,trial_function,a_in = float(a_range[0]),a_out = float(a_range[1]),k=1,a=1,mu=0)
        M = np.random.random()*np.pi*2
        
        sim.add(m=m,a=a,M=M,e=0,pomega=0,jacobi_masses=False)
    return sim

def GenerateEccentric(sim,N,a_range=[1,2],e_range=[0.0,0.1],massless=False):
    for i in tqdm(range(N)):
        if(massless):
            m = M_sun
            # r = 2*m/(1.37*10**4)**2 # c = 1.37*10^4 (radius units)/(time units) assuming
            # a period around central object is (2 * pi) = (1 year), (1 radius) = (1 au)
        else:
            m=0
            # r = 2*(1/332900)/(1.37*10**4)**2 # same as before but assuming mass of earth
        a = GenerateParameter(1,trial_function,a_in = float(a_range[0]),a_out = float(a_range[1]),k=1,a=1,mu=0)
        e = np.random.random()*(e_range[1]-e_range[0])+e_range[0]
        M = np.random.uniform(0,np.pi*2)
        pomega = np.random.uniform(0,np.pi*2)
        sim.add(m=m,a=a,M=M,e=e,pomega=pomega,jacobi_masses=False)
    return sim

# def GenerateEccentricFixed(sim,N):
#     for i in range(N):
#         m = 5e-5
#         a = 1.5
#         e = 0.3
#         inc = 0
#         Omega = np.random.uniform(0,np.pi*2)
#         M = np.random.uniform(0,np.pi*2)
#         omega = 0
#         sim.add(m=m,a=a,M=M,e=e,inc=inc,Omega=Omega,omega=omega,jacobi_masses=False)
#     return sim

def GenerateSphere(sim,N,massless=False):
    for i in tqdm(range(N)):
        if(massless):
            m = np.random.random()*1e-4
            # r = 2*m/(1.37*10**4)**2 # c = 1.37*10^4 (radius units)/(time units) assuming
            # a period around central object is (2 * pi) = (1 year), (1 radius) = (1 au)
        else:
            m=0
            # r = 2*(1/332900)/(1.37*10**4)**2 # same as before but assuming mass of earth
        a = GenerateParameter(1,trial_function,1,2,k=1/2,a=1,mu=0.999) # radius
        inc = np.random.random()*np.pi - np.pi/2 # inclincation above x-y plane, centered on x-axis
        M = np.random.uniform(0,np.pi*2) # true anomaly, where in the orbit the particle is
        Omega = np.random.uniform(0,np.pi*2) # angle from x-axis of ascending node
        sim.add(m=m,a=a,M=M,inc=inc,Omega=Omega,jacobi_masses=False)
    return sim

def GenerateEccentricSphere(sim,N,massless=False):
    for i in tqdm(range(N)):
        if(massless):
            m = np.random.random()*1e-4
            # r = 2*m/(1.37*10**4)**2 # c = 1.37*10^4 (radius units)/(time units) assuming
            # a period around central object is (2 * pi) = (1 year), (1 radius) = (1 au)
        else:
            m=0
            # r = 2*(1/332900)/(1.37*10**4)**2 # same as before but assuming mass of earth
        a = GenerateParameter(1,trial_function,1,2,k=1/2,a=1,mu=0.999) # radius
        e = np.random.random() # eccentricity, periapsis centered on x-axis
        inc = np.random.random()*np.pi - np.pi/2 # inclincation above x-y plane, centered on x-axis
        M = np.random.uniform(0,np.pi*2) # space ship position
        omega = np.random.uniform(0,np.pi*2) # angle from ascending node of periapsis
        Omega = np.random.uniform(0,np.pi*2) # angle from x-axis of ascending node
        sim.add(m=m,a=a,e=e,inc=inc,M=M,omega=omega,Omega=Omega,jacobi_masses=False)
    return sim

def GenerateCollision(sim,N,massless=False):
    for i in tqdm(range(N)):
        if(massless):
            m = M_sun
        else:
            m=0
        a = 1
        e = e_sun
        M = np.random.uniform(0,np.pi*2)
        pomega = np.random.uniform(0,np.pi*2)
        sim.add(m=m,a=a,e=e,M=M,pomega=pomega,r=R_sun,jacobi_masses=False)
    return sim

###########################################
# TOOLS
###########################################

def kick(sim,vx = 2/3):
    sim.move_to_hel()
    sim.particles[0].vx = vx
    sim.move_to_hel()
    return sim

def remove_param(sim,param="e",max_cut = False,minimum=1,maximum=100):
    length = len(sim.particles)
    reduc = 0
    dicty = {"e": [x.e for x in sim.particles[1:]],\
             "a": [x.a for x in sim.particles[1:]],\
             "inc": [x.inc for x in sim.particles[1:]],\
             "m": [x.m for x in sim.particles[1:]],\
             "omega": [x.omega for x in sim.particles[1:]],\
             "Omega": [x.Omega for x in sim.particles[1:]],\
             "M": [x.M for x in sim.particles[1:]],\
             "pomega": [x.pomega for x in sim.particles[1:]]}
        
    for i in range(1,length):
        if(max_cut == True):
            if(dicty[param][i-1] > maximum):
                sim.remove(i-reduc)
                reduc+=1
        else:
            if(dicty[param][i-1] < minimum):
                sim.remove(i-reduc)
                reduc+=1
                    
    return sim





####################################################
# NOW UNNESCECARRY
####################################################
def phys_shot(sim,filename):
    m_data = []
    x_data = []
    y_data = []
    z_data = []
    vx_data = []
    vy_data = []
    vz_data = []
    for i in range(len(sim.particles)):
        m_data.append(sim.particles[i].m)
        x_data.append(sim.particles[i].x)
        y_data.append(sim.particles[i].y)
        z_data.append(sim.particles[i].z)
        vx_data.append(sim.particles[i].vx)
        vy_data.append(sim.particles[i].vy)
        vz_data.append(sim.particles[i].vz)
    data = {'mass':m_data,'x':x_data,'y':y_data,'z':z_data,'vx':vx_data,'vy':vy_data,'vz':vz_data,'time':sim.t}
    df = pd.DataFrame(data)
    df.to_csv(filename)
    return sim

# also unnecessary
def unpack_phys_shot(sim,filename):
    df = pd.read_csv(filename)
    if(len(sim.particles) == 0):
        for i in range(df.shape[0]):
            sim.add(m=df['mass'][i],x=df['x'][i],y=df['y'][i],z=df['z'][i],vx=df['vx'][i],vy=df['vy'][i],vz=df['vz'][i])
        sim.t = df['time'][0]
    elif(df.size == len(sim.particles)):
        for i in range(len(sim.particles)):
            sim.particles[i].m = df['mass'][i]
            sim.particles[i].x = df['x'][i]
            sim.particles[i].y = df['y'][i]
            sim.particles[i].z = df['z'][i]
            sim.particles[i].vx = df['vx'][i]
            sim.particles[i].vy = df['vy'][i]
            sim.particles[i].vz = df['vz'][i]
        sim.t=df['time'][0]
    else:
        print("Sim is not empty and does not match the size of the imported sim data")
    return sim

def scrnsht_gen(sim,N,end,begin=0,dir=".",title="",Time = 0.5):
    sim.widget(size=(400,400))
    dest = dir+"/scrnshts_"+title
    timeseries = np.linspace(begin,end,N)
    for i in range(N):
        sim.output_screenshot(dest+f"/filename{i:0>5}.png")
        time.sleep(Time)
        sim.integrate(timeseries[i])
    # w.takeScreenshot(timeseries)
    # w.takeScreenshot(timeseries,prefix=dest+"/./screenshot"
    return sim

def simEvolve(sim,Time,steps,dt = 0.1):
    counter = 0
    for i in range(steps):
        sim.integrate(counter+dt)
        counter+=dt
        time.sleep(Time)
    print("Done")
    return sim