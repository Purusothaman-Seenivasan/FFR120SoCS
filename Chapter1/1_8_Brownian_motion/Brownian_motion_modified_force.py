import numpy as np
import random
from matplotlib import pyplot as plt

e0 = 0.0001
sigma0=1
mass0=1

velocity0 = np.sqrt((2*e0)/mass0)
t0= sigma0*np.sqrt(mass0/(2*e0))

deltaT=0.00005*t0
boundary=[0,100*sigma0]
n=50
iteration =int(50/deltaT)

int_velocity = np.full((n,2),20*velocity0)
print(int_velocity)
for i in range(n):
    angle= random.random()*2*np.pi
    int_velocity[i][0]= np.sin(angle)*int_velocity[i][0]
    int_velocity[i][1]=np.cos(angle)*int_velocity[i][1]
print(int_velocity)

LP_int_position= [boundary[1]/2, boundary[1]/2]
LP_int_velocity=np.array([0,0])
LP_radius= 10
LP_mass= 40*mass0

int_positions = []
def distance(point1, point2):
    return np.linalg.norm(point1 - point2)


def is_inside_circle(x, LP_int_position, LP_radius):
    return np.linalg.norm(x-LP_int_position)<=LP_radius

def distance_between(x, position):
    d=np.zeros((n,1))
    for i in range (n):
        d[i]= np.linalg.norm(x-position[i])-LP_radius
    return d

# def force(d):
#     f = np.zeros((n,1))
#     for i in range (n):
#         if d[i]>sigma0:
#             # f[i]= 4*e0*(((sigma0**12/d[i]**12)) - (sigma0**6/d[i]**6))
#             f[i]= 48*e0*((2*(sigma0**12/d[i]**13)) - (0.5*sigma0**6/d[i]**7))
#         else:
#             f[i]=0
#     return f

def lennard_jones_force(positions, x):
    forces = np.zeros_like(positions)
    for i in range(len(positions)):
        r = np.sqrt(np.sum((positions[i] - x) ** 2)) - LP_radius
        magnitude = 4 * e0 * (12 * np.power(sigma0, 12) * np.power(r, -13) - 6 * np.power(sigma0, 6) * np.power(r, -7))
        direction = (x - positions[i]) / r
        forces[i] = magnitude * direction

    return forces



while len(int_positions) < n:
    new_position = np.random.uniform(boundary[0], boundary[1], size=2)
    valid = True
    for existing_position in int_positions:
        if distance(new_position, existing_position) < sigma0 or is_inside_circle(new_position, LP_int_position, LP_radius):
            valid = False
            break
    if valid:
        int_positions.append(new_position)

int_positions= np.array(int_positions)
print(int_positions[:][0])
plt.xlim(0, boundary[1]) 
plt.ylim(0, boundary[1]) 
plt.quiver(int_positions[:, 0], int_positions[:, 1], int_velocity[:, 0], int_velocity[:, 1])

plt.show()       

for i in range(n):
    angle = np.random.uniform(0,2*np.pi) #random.random()*2*np.pi
    int_velocity[i][0]= int_velocity[i][0]*np.cos(angle) 
    int_velocity[i][1]= int_velocity[i][1]*np.sin(angle)


position= int_positions
velocity = int_velocity
Lp_position = LP_int_position
Lp_velocity = LP_int_velocity
Lp_position_list=[]
# Lp_position_list.append(Lp_position)
position_list=[]
# for t in range (iteration):
#     nxt_position = np.zeros((n,2))
#     nxt_velocity =np.zeros ((n,2))
    
#     Lp_next_position = np.zeros((1,2))
#     Lp_next_velocity = np.zeros((1,2))

#     Lp_middle= Lp_position + Lp_velocity*deltaT/2
#     d_LP= distance_between (Lp_middle, position) #calcualte distance between large particle and air
#     f_LP= lennard_jones_force(position, Lp_middle)
#     print(f_LP)
#     total_f=0#np.zeros((1,2))
#     for i in range (n):
#         total_f+= f_LP[i]

#     Lp_next_velocity = Lp_velocity+  total_f*deltaT/LP_mass
#     Lp_next_position = Lp_middle + Lp_next_velocity*deltaT/2

#     Lp_position= Lp_next_position
#     Lp_velocity = Lp_next_velocity
#     Lp_position_list.append(Lp_next_position)

for t in range(iteration):
    nxt_position = np.zeros((n, 2))
    nxt_velocity = np.zeros((n, 2))

    Lp_middle = Lp_position + 0.5 * Lp_velocity * deltaT
    f_LP = lennard_jones_force(position, Lp_middle)
    
    total_f = np.sum(f_LP, axis=0)  
    
    Lp_next_velocity = Lp_velocity + total_f * deltaT / LP_mass
    # Lp_next_position = Lp_middle+ Lp_next_velocity*deltaT/2
    Lp_next_position = Lp_position + Lp_velocity * deltaT + 0.5 * total_f * deltaT**2 / LP_mass

    if Lp_next_position[0]< boundary[0]+LP_radius:
        difference = boundary[0] - Lp_next_position[0]
        Lp_next_position[0]= difference
        Lp_next_velocity[0]= - Lp_next_velocity[0]
    if Lp_next_position[0]> boundary[1]-LP_radius:
        difference = Lp_next_position[0] - boundary[1]  
        Lp_next_position[0]= boundary[1] - difference
        Lp_next_velocity[0] = -Lp_next_velocity[0]
        
    if Lp_next_position[1]< boundary[0]+LP_radius:
        difference = boundary[0] - Lp_next_position[1]
        Lp_next_position[1]= difference
        Lp_next_velocity[1]= - Lp_next_velocity[1]
    if Lp_next_position[1]> boundary[1]-LP_radius:
        difference = Lp_next_position[1] - boundary[1]  
        Lp_next_position[1]= boundary[1] - difference
        Lp_next_velocity[1] = -Lp_next_velocity[1] 
    Lp_position = Lp_next_position
    Lp_velocity = Lp_next_velocity
    Lp_position_list.append(Lp_next_position)

    for i in range (n):
        nxt_position[i] = position[i]+ velocity[i]*deltaT
        nxt_velocity[i] = velocity[i]

        if nxt_position[i][0] < boundary[0] or nxt_position[i][0] > boundary[1]:
            nxt_position[i][0] = max(min(nxt_position[i][0], boundary[1]), boundary[0])
            nxt_velocity[i][0] *= -1  # Reverse the velocity for reflection

        # Boundary reflection for y-coordinate
        if nxt_position[i][1] < boundary[0] or nxt_position[i][1] > boundary[1]:
            nxt_position[i][1] = max(min(nxt_position[i][1], boundary[1]), boundary[0])
            nxt_velocity[i][1] *= -1
    position=nxt_position
    velocity = nxt_velocity
    position_list.append(position)       
        # if is_inside_circle (nxt_position[i], Lp_position, LP_radius):
        #     direction = nxt_position[i] -Lp_position
        #     normalized_direction = direction/np.linalg.norm(direction)
        #     nxt_position[i]= Lp_position+normalized_direction*LP_radius
        #     nxt_velocity[i] = nxt_velocity[i]-2*np.dot(nxt_velocity[i],normalized_direction)*normalized_direction






    



# circle = plt.Circle(LP_int_position, LP_radius, edgecolor='black', facecolor='none')
# plt.figure()
# plt.gca().set_aspect('equal', adjustable='box')
# # plt.xlim(boundary)
# # plt.ylim(boundary)
# plt.gca().add_patch(circle)
# plt.quiver(LP_int_position[0], LP_int_position[1], LP_int_velocity[0], LP_int_velocity[1])

# for i in range(n):
#     plt.scatter(int_positions[i][0], int_positions[i][1], color='blue')    
# plt.xlabel('X-axis')
# plt.ylabel('Y-axis')
# plt.title('Circle and Particles Plot')
# plt.grid(True)
# plt.show()
# #------------------------------
position_list= np.array(position_list)
Lp_position_list=np.array(Lp_position_list)


# circle = plt.Circle(Lp_position, LP_radius, edgecolor='black', facecolor='none')
# plt.figure()
# plt.gca().set_aspect('equal', adjustable='box')
# plt.gca().add_patch(circle)
# plt.scatter(position[:, 0], position[:, 1])
# # plt.quiver(position[:, 0], position[:, 1], velocity[:, 0], velocity[:, 1], color='r')

# # plt.plot(Lp_position_list[:,0],Lp_position_list[:,0], label=f"Particle")

# plt.show()

#----------------------------------------



# Plotting the circle and particles
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

# Plotting Circle and Particles
circle = plt.Circle(LP_int_position, LP_radius, edgecolor='black', facecolor='none')
ax1.set_aspect('equal')
ax1.add_patch(circle)
ax1.quiver(LP_int_position[0], LP_int_position[1], LP_int_velocity[0], LP_int_velocity[1])
for i in range(n):
    ax1.scatter(int_positions[i][0], int_positions[i][1], color='blue')

ax1.set_xlabel('X-axis')
ax1.set_ylabel('Y-axis')
ax1.set_title('Circle and Particles Plot')
ax1.grid(True)

# Plotting the trace of a particle
ax2.set_aspect('equal')
circle = plt.Circle(Lp_position, LP_radius, edgecolor='black', facecolor='none')
ax2.add_patch(circle)
ax2.scatter(position[:, 0], position[:, 1])
ax2.plot(Lp_position_list[:, 0], Lp_position_list[:, 1], label='Trace of Particle')
ax2.set_xlabel('X-axis')
ax2.set_ylabel('Y-axis')
ax2.set_title('Trace of Particle')
ax2.grid(True)

plt.tight_layout()
plt.show()


msd= np.zeros((iteration,1))
iteration_list=[]
for t in range (iteration):
    dif_square=0
    for i in range(iteration-t):
        dif_square+= np.square(np.linalg.norm(Lp_position_list[i+t]-Lp_position_list[i]))
    msd[t]= 1/(iteration-t) * dif_square
    iteration_list.append(t)

slope, _ = np.polyfit(iteration_list, msd, 1)
# conversion_factor = 1e-6 
diffusion_coefficient = slope / 4
print('diff', diffusion_coefficient)
# D = msd/(4*1)
# print(D)
plt.plot(iteration_list,msd, 'o')
plt.plot(iteration_list, slope*iteration_list, '-', color='black' )
plt.xlim(0,iteration)
plt.show()

