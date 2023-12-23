from matplotlib import pyplot as plt
import numpy as np

intial_position = 0.0
intial_velocity = 0.1
time=0
delta_t = 0.002
spring_constant = 2 #k
mass = 1

time_list=[]
time_list.append(time)
iteration = int(10/delta_t)

positions = np.zeros(iteration)
velocity = np.zeros(iteration)
positions[0]= intial_position
velocity[0] = intial_velocity

total_energy = np.zeros(iteration)
# # potential_energy = 1/2*spring_constant*position**2
# kinetic_energy = 1/2*mass* velocity**2
# total_energy = potential_energy + kinetic_energy
# total_energy_list=[]
# total_energy_list.append(total_energy)


for i in range (iteration-1):
    omega_square = spring_constant/ mass 
    a= -omega_square* positions[i]
    positions[i+1] = positions[i] + velocity[i] * delta_t
    # velocity = velocity + a*delta_t
    velocity[i+1] = velocity[i] + a*delta_t
    
    time = time + delta_t
    time_list.append(time)

for i in range(iteration):
    kinetic_energy= 0.5*mass*velocity[i]**2
    potential_energy=  0.5*spring_constant*positions[i]**2
    total_energy[i]= kinetic_energy+potential_energy

    # kinetic_energy[i+1] = 0.5*mass* velocity[i]**2
    # potential_energy = 0.5*spring_constant*position[i]**2
    # total_energy = potential_energy + kinetic_energy 
    # total_energy_list.append(total_energy)
# for i in range(iteration):
#----------Leap frog method--------------

lf_position=np.zeros(iteration)
lf_velocity = np.zeros(iteration)
lf_total_energy =np.zeros(iteration)

lf_position[0]= intial_position
lf_velocity[0]= intial_velocity
print('intial', lf_position)
for i in range (iteration-1):
    middle_step= lf_position[i]+lf_velocity[i]*delta_t/2
    print('middle', middle_step)
    a= - (spring_constant/mass) *middle_step
    lf_velocity[i+1]= lf_velocity[i]+ a*delta_t
    lf_position[i+1] = middle_step  + lf_velocity[i+1]*delta_t/2


for i in range(iteration):
    lf_kinetic_energy= 0.5*mass*lf_velocity[i]**2
    lf_potential_energy=  0.5*spring_constant*lf_position[i]**2
    lf_total_energy[i]= lf_kinetic_energy+lf_potential_energy

continous_time = np.arange(0, 10, 0.01)
print(len(continous_time))

r0= 0.1
v0=0
omega_alalytical= np.sqrt(spring_constant/mass)
a_analytical = np.sqrt(intial_position**2+(intial_velocity/omega_alalytical)**2)
phi = np.arctan2((-intial_velocity/(a_analytical*omega_alalytical)),(intial_position/a_analytical))


total_energy_analytical_value =[]
functionvalues=[]
velocity_analytical_list=[]

for t in time_list:
    position_function = a_analytical*(np.cos(omega_alalytical*t+phi))
    functionvalues.append(position_function)

    velocity_function = -omega_alalytical*a_analytical*np.sin(omega_alalytical*t+phi)
    velocity_analytical_list.append(velocity_function)
    
    potential_energy_analytical= 0.5*mass*omega_alalytical**2*a_analytical**2*np.sin(omega_alalytical*t+phi)**2
    kinetic_energy_analytical = 0.5*spring_constant*a_analytical**2*np.cos(omega_alalytical*t+phi)**2
    total_energy_analytical = potential_energy_analytical+ kinetic_energy_analytical
    # total_energy_analytical = 0.5*spring_constant*a_analytical**2
    total_energy_analytical_value.append(total_energy_analytical)

    # potential_energy = 1/2*spring_constant*r0**2
    # kinetic_energy = 1/2*mass* (omega*np.sqrt(a**2-))**2
    # total_energy = potential_energy + kinetic_energy


plt.plot(time_list,functionvalues , color='r', linestyle ='-' ,label='analytical solution of Euler algorithm- position' )
plt.plot(time_list,positions, color='g', linestyle= '--' ,label='numerical solution of Euler algorithm - position')
plt.legend()
plt.show()  

plt.plot(time_list,velocity_analytical_list , color='r', linestyle ='-' ,label='analytical solution of Euler algorithm- velocity' )
plt.plot(time_list,velocity, color='g', linestyle= '--' ,label='numerical solution of Euler algorithm - velocity')
plt.legend()
plt.show() 

plt.plot(time_list,total_energy, color='r', linestyle ='-',label='analytical solution of Euler algorithm - Total Energy')
plt.plot(time_list,total_energy_analytical_value, color='g', linestyle= '--' ,label='numerical solution of Euler algorithm - Total Energy')
plt.legend()
plt.show()


plt.plot(time_list,positions, color='b', linestyle= 'dotted' ,label='numerical solution of Euler algorithm - position')
plt.plot(time_list,functionvalues , color='r', linestyle ='-' ,label='analytical solution- position' )
plt.plot(time_list,lf_position, color='g', linestyle= '--' ,label='numerical solution of leaf frog algorithm - position')
plt.legend()
plt.show()  


plt.plot(time_list,total_energy, color='b', linestyle ='dotted',label='analytical solution of Euler algorithm - Total Energy')
plt.plot(time_list,lf_total_energy, color='r', linestyle ='-',label='analytical solution of Leaf frog algorithm - Total Energy')
plt.plot(time_list,total_energy_analytical_value, color='g', linestyle= '--' ,label='numerical solution of Euler algorithm - Total Energy')
plt.legend()
plt.show()

# print(iteration)
# print(time_list)