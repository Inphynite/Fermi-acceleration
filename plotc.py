import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import math

# This Python script plots data read on a txt file, generated from code written in a C script

df = pd.read_csv("data.txt")


u = df["vit"]
psi = df["ph"]


# Conversion from pandas DataFrames to numpy arrays (converting string to float)
u2 = np.array(u)
u3 = u2.astype(float)
psi2 = np.array(psi)
psi3 = psi2.astype(float)


print(max(u3))
# Plot
plt.figure(figsize=(5, 8))

plt.plot(psi3, u3, ',')

# For zooming in on the maps
#plt.ylim(8.75, 9)
#plt.xlim(2.9, 3.3)


plt.show()

# Histogram for velocity distribution
#plt.hist(u3, bins=200)
plt.hist(u3, bins=100, log=True)
plt.xlabel("u")
plt.ylabel("Nombre d'occurrences")
#plt.xscale("log") # For log scale on both axes, doesn't look quite right, different size bins
plt.show()

# Velocity distribution using seaborn
#sb.kdeplot(u3, bw_adjust=0.3)


# Plot of boundary velocities as a function of M (boundary velocities have been drawn from the same C script mentioned in the beginning, run with different values of M)
M = np.array([10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
umax = np.array([12.806066, 12.318397, 16.800583, 22.408844, 23.622002, 27.508723, 32.444714, 30.756882, 34.538208, 33.83948, 36.109803])**2
x = np.linspace(10, 200, 50)
y = 2.65 * np.sqrt(x)
plt.plot(x, y**2, label='Ajustement')
plt.plot(M, umax, '.', label='Valeurs calculees')
plt.xlabel("M")
plt.ylabel("$u_{b}^{2}$")
plt.legend()
plt.show()
