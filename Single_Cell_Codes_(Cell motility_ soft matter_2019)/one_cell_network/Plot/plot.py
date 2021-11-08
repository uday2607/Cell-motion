import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from math import pi, cos, sin
plt.style.use('seaborn-pastel')
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
FFwriter = animation.FFMpegWriter(fps = 1)


t = np.linspace(0, 2*pi, 100)

Lx = 60
Ly = 60
a0 = 4
b = 4

x = []
y = []
with open("testf10id0xc.txt", "r") as f:
    lines = f.readlines()

for l in lines:
    x_, y_ = [abs(float(i)%Lx) for i in l.split(" ")]
    x.append(x_)
    y.append(y_)

theta = []
with open("testf10id0theta.txt", "r") as f:
    lines = f.readlines()

for l in lines:
    t_ = float(l)
    theta.append(t_*pi/180)

print(x,y)
fig = plt.figure()
ax = plt.axes(xlim=(0, Lx), ylim=(0, Ly))
line, = ax.plot([], [], lw=2)

# initialization function
def init():
	# creating an empty plot/frame
    a = a0
    Ell = np.array([a*np.cos(t) , b*np.sin(t)])
    R_rot = np.array([[cos(theta[0]) , -sin(theta[0])],[sin(theta[0]) , cos(theta[0])]])
    Ell_rot = np.zeros((2,Ell.shape[1]))
    for i in range(Ell.shape[1]):
        Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])
    line.set_data(x[0]+Ell_rot[0,:] , y[0]+Ell_rot[1,:])
    return line,

# animation function
def animate(i):

    a = a0*(1-0.4*(i%30)/30)
    Ell = np.array([a*np.cos(t) , b*np.sin(t)])
    R_rot = np.array([[cos(theta[i]) , -sin(theta[i])],[sin(theta[i]) , cos(theta[i])]])
    Ell_rot = np.zeros((2,Ell.shape[1]))
    for j in range(Ell.shape[1]):
        Ell_rot[:,j] = np.dot(R_rot,Ell[:,j])
    line.set_data(x[i]+Ell_rot[0,:] , y[i]+Ell_rot[1,:])
    return line,

# setting a title for the plot
plt.title('Cell Motion')
# hiding the axis details
plt.axis('off')
plt.grid(color='lightgray',linestyle='--')

# call the animator
anim = animation.FuncAnimation(fig, animate, init_func=init,
							frames=len(x), interval=200, blit=True)
anim.save('animation.mp4', writer = FFwriter)

# save the animation as mp4 video file
plt.show()
