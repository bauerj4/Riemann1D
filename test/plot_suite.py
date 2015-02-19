import numpy as np
import os
import sys
import matplotlib.pyplot as plt

PATH = sys.argv[1];

f = open(PATH);

position = [];
velocity = [];
density = [];
pressure = [];


line = f.readline();

while(line != ""):
    tokens = line.split();
    #for x in tokens:
    #    print x
    position.append(float(tokens[0]));
    density.append(float(tokens[1]));
    pressure.append(float(tokens[2]));
    velocity.append(float(tokens[3]));
    line = f.readline();

postion = np.array(position);
velocity = np.array(velocity);
pressure = np.array(pressure);
density = np.array(density);

plt.scatter(position,velocity);
plt.xlabel("Position");
plt.ylabel("Velocity");
plt.xlim(np.amin(position), np.amax(position))
#plt.ylim(1.25 * np.amin(velocity), 1.25 * np.amax(velocity));
plt.show();

plt.scatter(position, pressure);
plt.xlabel("Position");
plt.ylabel("Pressure");
plt.xlim(np.amin(position), np.amax(position))
plt.ylim(0., 1.25 * np.amax(pressure));
plt.show();

plt.scatter(position, density);
plt.xlabel("Position");
plt.ylabel("Density");
plt.xlim(np.amin(position), np.amax(position))
plt.ylim(0., 1.25 * np.amax(density));
plt.show();
