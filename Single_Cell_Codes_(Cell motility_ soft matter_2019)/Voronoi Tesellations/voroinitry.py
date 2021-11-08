from numpy import *
from pylab import *
from scipy.spatial import Voronoi, voronoi_plot_2d
from random import random
points=[]
N=400; RN=sqrt(N)
for i in arange(RN):
    for j in arange(RN):
        x=i+0.5*(random()-0.5)
        y=j+0.5*(random()-0.5)
        r1=2*pi*random()
        dx=0.25*cos(r1); dy=0.25*sin(r1)
        points=points+[[x+dx,y+dy],[x-dx,y-dy]]
points=asarray(points)
plot(points[:,0],points[:,1],'w.')
axis('square')
vor=Voronoi(points)
#voronoi_plot_2d(vor)
skmax=shape(vor.ridge_vertices)[0]
for sk in arange(skmax):
    #dindex=abs(vor.ridge_points[sk,0]-vor.ridge_points[sk,1])
    a=vor.ridge_points[sk,0]
    b=vor.ridge_points[sk,1]
    c=min(a,b); d=abs(a-b)
    flag=1
    if d==1 and mod(c,2)==0: flag =0
    if flag==1:
        s=vor.ridge_vertices[sk]
        s= asarray(s)
        if all(s>=0):
            plot(vor.vertices[s,0],vor.vertices[s,1],'k-')
for k in arange(0,2*N,2):
    plot([points[k,0],points[k+1,0]],[points[k,1],points[k+1,1]],'r-')
    
show()

        
