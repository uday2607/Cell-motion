import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
plt.style.use('seaborn-pastel')

x = []
y = []
with open("xc.txt", "r") as f:
    lines = f.readlines()

for l in lines:
    x_, y_ = [float(i) for i in l.split(" ")]
    x.append(x_)
    y.append(y_)

fig = plt.figure()
ax = plt.axes(xlim=(35, 36), ylim=(25, 26))
line, = ax.plot([], [], lw=2)

# initialization function
def init():
	# creating an empty plot/frame
	line.set_data([x[1], x[0]], [y[1], y[0]])
	return line,

# animation function
def animate(i):

    line.set_data([x[i], x[i-1]], [y[i], y[i-1]])
    return line,

# setting a title for the plot
plt.title('Cell Motion')
# hiding the axis details
plt.axis('off')

# call the animator
anim = animation.FuncAnimation(fig, animate, init_func=init,
							frames=len(x), interval=200, blit=True)

# save the animation as mp4 video file
plt.show()
