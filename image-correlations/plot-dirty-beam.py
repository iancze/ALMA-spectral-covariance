import numpy as np 
from astropy.io import ascii
import matplotlib.pyplot as plt 

from dataset import vis_correlated as vis

arcsec = np.pi / (180.0 * 3600)

# # augment Hermitian pairs.
# V = d.V_reals + 1.0j * d.V_imags

# us = np.concatenate([d.u, -d.u]) * 1e3 # shape (2 * nsamples)
# vs = np.concatenate([d.v, -d.v]) * 1e3 # shape (2 * nsamples)
# Vs = np.concatenate([V, np.conj(V)], axis=1) # shape (nchan, 2 * nsamples)
# ws = d.weight * np.ones_like(us) # shape (2 * nsamples)

# w_2 = d.weight * np.ones_like(d.u)



# calculate the normalization constant
C = 1 / np.sum(vis.weightweight)

dl = 0.02 # arcseconds
npix = 256
ls = dl * (np.arange(npix) - npix/2)
ms = dl * (np.arange(npix) - npix/2)

# create a 2D grid 
ll, mm = np.meshgrid(ls * arcsec, ms * arcsec) # radians
# shape (npix, npix)

# Choose a channel to image
chan = vis.nchan//2

# calculate beam values over grid 
# broadcast loose visibility arrays to each and every l,m cell.
beam = C * np.sum(vis.weightweight * np.exp(2.0j * np.pi * (1e3 * vis.uu[np.newaxis, np.newaxis, :] * ll[:, :, np.newaxis] + 1e3 * vis.vv[np.newaxis, np.newaxis, :] * mm[:, :, np.newaxis])), axis=2)

# (left, right, bottom, top),
extent = [np.max(ls) + 0.5 * dl, np.min(ls) - 0.5 * dl, np.min(ms) - 0.5 * dl, np.max(ms) + 0.5 * dl]

fig, ax = plt.subplots(nrows=1, figsize=(5,5))
im = ax.imshow(np.fliplr(beam.real), origin="lower", extent=extent, cmap="inferno")
ax_cb = plt.colorbar(im, ax=ax)
ax_cb.set_label(r"$B_D$")
ax.set_xlabel(r"$\Delta \alpha \cos \delta$ [arcseconds]")
ax.set_ylabel(r"$\Delta \delta$ [arcseconds]")
fig.subplots_adjust(top=0.98, left=0.2)

fig.savefig("beam.png", dpi=200)


# calculate our spatial covariance function 
# set l_1 = 0, m_1 = 0
spatial_cov = 4 * C**2 * np.sum(vis.weight * np.cos(2 * np.pi * (1e3 * vis.u[np.newaxis, np.newaxis, :] * ll[:, :, np.newaxis] + 1e3 * vis.v[np.newaxis, np.newaxis, :] * mm[:, :, np.newaxis])), axis=2) 

fig, ax = plt.subplots(nrows=1, figsize=(5,5))
im = ax.imshow(np.fliplr(spatial_cov), origin="lower", extent=extent, cmap="inferno")
ax_cb = plt.colorbar(im, ax=ax)
ax_cb.set_label(r"$k$")
ax.set_xlabel(r"$\Delta \alpha \cos \delta$ [arcseconds]")
ax.set_ylabel(r"$\Delta \delta$ [arcseconds]")
fig.subplots_adjust(top=0.98, left=0.2)

fig.savefig("spatial-cov.png", dpi=200)


# I think the spatial covariance function and the dirty beam end up being the same as each other
# because the dirty beam autocorrelated with itself is still the dirty beam. That means that even
# though the spatial covariance function has units of ^2, it has the same shape as the dirty beam.
fig, ax = plt.subplots(nrows=1, figsize=(5,5))
im = ax.imshow(np.fliplr(beam.real) - np.fliplr(spatial_cov)/np.max(spatial_cov), origin="lower", extent=extent, cmap="inferno")
ax_cb = plt.colorbar(im, ax=ax)
ax_cb.set_label(r"$k$")
ax.set_xlabel(r"$\Delta \alpha \cos \delta$ [arcseconds]")
ax.set_ylabel(r"$\Delta \delta$ [arcseconds]")
fig.subplots_adjust(top=0.98, left=0.2)

fig.savefig("difference.png", dpi=200)
