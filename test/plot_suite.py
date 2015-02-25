import numpy as np
import os
import sys
import matplotlib.pyplot as plt


GAMMA = 1.4;
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

internal_energy = pressure/( (GAMMA -1) * density);

plt.figure(1)

plt.subplot(221)
plt.scatter(position,velocity,color='black',marker='.');
plt.xlabel("Position");
plt.ylabel("Velocity");
plt.xlim(np.amin(position), np.amax(position))
#plt.ylim(1.25 * np.amin(velocity), 1.25 * np.amax(velocity));
#plt.show();

plt.subplot(222)

plt.scatter(position, pressure,color='black',marker='.');
plt.xlabel("Position");
plt.ylabel("Pressure");
plt.xlim(np.amin(position), np.amax(position))
plt.ylim(0., 1.25 * np.amax(pressure));
#plt.show();

plt.subplot(223)

plt.scatter(position, density,color='black',marker='.');
plt.xlabel("Position");
plt.ylabel("Density");
plt.xlim(np.amin(position), np.amax(position))
plt.ylim(0., 1.25 * np.amax(density));
#plt.show();

plt.subplot(224)

plt.scatter(position, internal_energy,color='black',marker='.');
plt.xlabel("Position");
plt.ylabel("Internal Energy");
plt.xlim(np.amin(position), np.amax(position))
plt.ylim(0., 1.25 * np.amax(internal_energy));
plt.show();
