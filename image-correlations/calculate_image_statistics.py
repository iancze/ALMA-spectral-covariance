import numpy as np 
import matplotlib.pyplot as plt 

import dataset as d

arcsec = np.pi / (180.0 * 3600)

# augment Hermitian pairs.
V = d.V_reals + 1.0j * d.V_imags

us = np.concatenate([d.u, -d.u]) * 1e3 # shape (2 * nsamples)
vs = np.concatenate([d.v, -d.v]) * 1e3 # shape (2 * nsamples)
Vs = np.concatenate([V, np.conj(V)], axis=1) # shape (nchan, 2 * nsamples)
ws = d.weight * np.ones_like(us) # shape (2 * nsamples)

# calculate the normalization constant
C = 1 / np.sum(ws)
w_2 = d.weight * np.ones_like(d.u)

# load the image cube we already calculated
dcube = np.load("image_cube.npz")
img_cube = dcube["img_cube"]
ls = dcube["ls"] # arcsec
ms = dcube["ms"] # arcsec
nchan = len(img_cube)

# calculate the variance of each pixel in the image, in just one channel
var = np.var(img_cube[0])
print("empricial var", var)

# calculate this theoretically from our expression
var_vis = 4 * C**2 * np.sum(w_2) # per-channel
print("theoretical var", var_vis)

print(var_vis / var)

# 11/27/22: We aren't able to get these expressions to match very well, we seem to be off by factors of 2 - 6 in the variance, depending on how we re-arrange sums and pre-factors.
# 11/28/22: We realized that the Hann-smoothing operation had changed the noise level, so we needed to rescale the weights.


# # calculate the pixel-to-pixel covariance within an image.
# # we can take a cheap way out and make this a 1D slice plot, assuming it will be similar in all directions.
# cov = np.cov(img_cube[0], rowvar=True)

# fig, ax = plt.subplots(nrows=1)
# ax.imshow(cov)
# fig.savefig("cov-matrix.png", dpi=200)

