import time
start_time = time.time()
import Functions
import numpy as np
import sys
import glob

############################
#input parameters are as follows
#------------------------------
# no inputs, default values specified in .py
# input 1, output csv filename
# input 2, input csv filename for starting simulation conditions
# input 3, kick
# input 4, time to integrate until
# input 5, number of particles
# input 6, mode
# input 7, massless before kick

# out1.bin None True 100 100 disk True

# valid modes: 'disk', 'eccen', 'sphere', 'eccen_sphere'
############################
# currently require one of the following
#  - 1 inputs
#  - 2 inputs
#  - 3 inputs
#  - 4 inputs
#  - 5 inputs
#  - 6 inputs
#  - 7 inputs

output = 'out' # 1st input
iinput = 'None' # 2nd input
kick = False # 3rd input
t = 10 # 4th input
N=1 # 5th input
mode = 'disk' # 6th input
massless = False # 7th input

#global variables
inputted = False
arguements = []

########Checking inputs###########
# no inputs
if(len(sys.argv)<2):
    print("Using default values")

#1 input: output
if(len(sys.argv)>=2):
    output = sys.argv[1]
    arguements.append(output)

#2 inputs: output input
if(len(sys.argv)>=3):
    iinput = sys.argv[2]
    arguements.append(iinput)

#3 inputs: output input kick
if(len(sys.argv)>=4):
    if(sys.argv[3] == 'True'):
        kick=True
    else:
        kick=False
    arguements.append(kick)

#4 inputs: output input kick time
if(len(sys.argv)>=5):
    t = float(sys.argv[4])
    arguements.append(t)

#5 inputs: output input kick time #_of_particles
if(len(sys.argv)>=6):
    N = int(sys.argv[5])
    arguements.append(N)

#6 inputs: output input kick time #_of_particles generation_method
if(len(sys.argv)>=7):
    mode = sys.argv[6]
    arguements.append(mode)

if(len(sys.argv)==8):
    if(sys.argv[7] == 'True'):
        massless=True
    else:
        massless=False
    arguements.append(massless)

# too many
if(len(sys.argv)>7):
    print('too many arguments so far, ignoring extras')

print(f"output filename = {output}, input filename = {iinput}, is kicking true?: {kick}\nintegrating time = {t}, number of particles = {N}, mode = {mode}\
, massless before kick?: {massless}")

#########Initialize simulation###########
mass_array = []
if(iinput == 'None'):
    print("Generating Particles")
    sim = Functions.Initialize()
    if(mode == 'disk'):
        sim = Functions.GenerateDisk(sim,N,massless=massless)
    elif(mode == 'eccen'):
        sim = Functions.GenerateEccentric(sim,N,massless=massless)
    elif(mode == 'sphere'):
        sim = Functions.GenerateSphere(sim,N,massless=massless)
    elif(mode == 'eccen_sphere'):
        sim = Functions.GenerateEccentricSphere(sim,N,massless=massless)
    else:
        sys.exit("Not a valid mode, try again")
else:
    inputted = True
    print("inputting particles")
    if(iinput[-4:]!='.bin'):
        sys.exit('not a .bin file for input, try again')
    # sim = Functions.Initialize(massed=False)
    # sim = Functions.unpack_phys_shot(sim,input)
    sim = rebound.Simulation(iinput)
    if(massless):
        for i in range(len(sim.particles)):
            mass_array.append(sim.particles[i].m)
            sim.particles[i].m = 0

sim.save_to_file(output[:-4]+'bk.bin')
# sim = Functions.phys_shot(sim,output[:-4]+'bk.csv')
###########Check kick############
if(kick):
    print("kicking")
    sim.move_to_hel()
    sim.particles[0].vx=1/2
    sim.move_to_hel()
    if(massless): # adding mass, this may be moved later
        if(inputted):
            for i in range(len(mass_array)):
                sim.particles[i].m = mass_array[i]
        else:
            for i in range(len(sim.particles)):
                m = np.random.random()*1e-4
                sim.particles[i].m = m

sim.save_to_file(output[:-4]+'dk.bin')
# sim = Functions.phys_shot(sim,output[:-4]+'dk.csv')
###########Integrate############
print("integrating")
sim.integrate(t)

###########Output csv###########
print("outputting")
sim.save_to_file(output[:-4]+'ak.bin')
# sim = Functions.phys_shot(sim,output[:-4]+'ak.csv')

print("--- %s seconds ---" % (time.time() - start_time))