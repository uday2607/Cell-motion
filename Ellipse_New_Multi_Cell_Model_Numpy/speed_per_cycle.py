import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

with open("distance_data.txt", "r") as f:
    lines = f.readlines()

t1 = 0
t2 = 0
dis = 0
x2, x1 = 0, 0
y2, y1 = 0, 0
distances = []
displacements = []
for line in lines[1: ]:
    x, y, t, p = line.strip("\n").split("\t")
    x = float(x)
    y = float(y)
    t = float(t)

    if t == 0 and p == "c":
        t1 = 0
        x1 = x
        y1 = y

    if p == "n":
        t2 = t 
        distances.append(dis/(t2 - t1))
        t1 = t2 
        dis = 0

    x2 = x
    y2 = y 
    dis += np.sqrt((x2-x1)**2.0 + (y2-y1)**2.0)    
    x1 = x2 
    y1 = y2    

for line in lines[1: ]:
    x, y, t, p = line.strip("\n").split("\t")
    x = float(x)
    y = float(y)
    t = float(t)

    if t == 0 and p == "c":
        t1 = 0
        x1 = x
        y1 = y

    if p == "n":
        x2 = x
        y2 = y 
        t2 = t 
        displacements.append(np.sqrt((x2-x1)**2.0 + (y2-y1)**2.0)/(t2 - t1))
        x1 = x2 
        y1 = y2 
        t1 = t2 
        dis = 0
    

distances = np.array(distances)
q25, q75 = np.percentile(distances, [25, 75])
bin_width = 2 * (q75 - q25) * len(distances) ** (-1/3)
bins = round((distances.max() - distances.min()) / bin_width)
print("Freedman–Diaconis number of bins:", bins)
sns.histplot(distances, kde=True, 
             bins=bins, color = 'darkblue')
plt.ylabel('Density')
plt.xlabel('Speed per cycle (bl/unit time)')
plt.show()

displacements = np.array(displacements)
q25, q75 = np.percentile(displacements, [25, 75])
bin_width = 2 * (q75 - q25) * len(displacements) ** (-1/3)
bins = round((displacements.max() - displacements.min()) / bin_width)
print("Freedman–Diaconis number of bins:", bins)
sns.histplot(displacements, kde=True, 
             bins=bins, color = 'darkblue')
plt.ylabel('Density')
plt.xlabel('Velocity per cycle (bl/unit time)')
plt.show()