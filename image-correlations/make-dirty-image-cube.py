import numpy as np 
from dataset import vis_correlated, vis_uncorrelated

arcsec = np.pi / (180.0 * 3600)

# # augment Hermitian pairs.
# V = d.V_reals + 1.0j * d.V_imags

# us = np.concatenate([d.u, -d.u]) * 1e3 # shape (2 * nsamples)
# vs = np.concatenate([d.v, -d.v]) * 1e3 # shape (2 * nsamples)
# Vs = np.concatenate([V, np.conj(V)], axis=1) # shape (nchan, 2 * nsamples)
# ws = d.weight * np.ones_like(us) # shape (2 * nsamples)


def create_cube(fname, vis):
    # calculate the normalization constant, per-channel
    C = 1 / np.sum(vis.weightweight)

    # set image parameters
    dl = 0.025 # arcseconds
    npix = 256
    ls = dl * (np.arange(npix) - npix/2) # [arcsec]
    ms = dl * (np.arange(npix) - npix/2)

    # create a 2D grid 
    ll, mm = np.meshgrid(ls * arcsec, ms * arcsec) # [radians]
    # shape (npix, npix)

    img_cube = np.empty((vis.nchan, *ll.shape))

    # calculate beam values over grid 
    # broadcast loose visibility arrays to each and every l,m cell.
    for i in range(vis.nchan):
        img_cube[i] = (C * np.sum(vis.weightweight * vis.VV[i, :] * np.exp(2.0j * np.pi * (1e3*vis.uu[np.newaxis, np.newaxis, :] * ll[:, :, np.newaxis] + 1e3*vis.vv[np.newaxis, np.newaxis, :] * mm[:, :, np.newaxis])), axis=2)).real

    np.savez(fname, img_cube=img_cube, ls=ls, ms=ms)

create_cube("image_cube_uncorrelated.npz", vis_uncorrelated)
create_cube("image_cube_correlated.npz", vis_correlated)