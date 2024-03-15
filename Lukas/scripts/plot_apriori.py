import numpy as np
import pickle
import sys
from matplotlib import pyplot as plt

profiles = []
for fname in sys.argv[1:]:
    with open(fname, "rb") as file:
        profiles.append(pickle.load(file))

plt.figure()
for i in range(len(profiles)):
    plt.plot(profiles[i][0], profiles[i][1] * 0.001, linestyle='-', label=sys.argv[i + 1])
    plt.plot(profiles[i][2], profiles[i][1] * 0.001, linestyle=':', label=sys.argv[i + 1])

plt.xlabel("VER")
plt.ylabel("Altitude, km")
plt.legend()
plt.show()
