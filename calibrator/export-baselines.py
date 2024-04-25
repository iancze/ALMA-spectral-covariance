import numpy as np
import matplotlib.pyplot as plt 

import casatools

msmd = casatools.msmetadata()
ms = casatools.ms()

fname = "calibrator_data_spw25_FDM.ms"
# msmd.open(fname)
# spws = msmd.datadescids() 
# print(spws)
# msmd.done()

spw_id = 0


ms.open(fname)
ms.selectinit(spw_id)
d = ms.getdata(["uvw"], ifraxis=True)  
ms.done()

# print(d["uvw"].shape)
# # d["uvw"] is an array of float64 with shape [3, nvis, ntime]
uu, vv, ww = d["uvw"]  # unpack into len nvis vectors

u = uu[:,0]
v = vv[:,0]

fig, ax = plt.subplots(nrows=1, figsize=(3.5, 3.5))
ax.scatter(u, v, s=1.5, rasterized=True, linewidths=0.0, c="k")
ax.set_xlabel(r"$u$ [m]")
ax.set_ylabel(r"$v$ [m]")

fig.savefig("baselines_m.png", dpi=300)

np.savez("snapshot-baselines.npz", u=u, v=v) # [m]
