import matplotlib.pyplot as plt 
import numpy as np

d = np.load("snapshot-baselines.npz")
u = d["u"]
v = d["v"]

fig, ax = plt.subplots(nrows=1, figsize=(3.5, 3.5))
ax.scatter(u, v, s=1.5, rasterized=True, linewidths=0.0, c="k")
ax.set_xlabel(r"$u$ [m]")
ax.set_ylabel(r"$v$ [m]")

fig.savefig("baselines_m.png", dpi=300)

