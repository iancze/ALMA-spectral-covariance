import numpy as np
from astropy.constants import c
from dataclasses import dataclass

@dataclass
class Visibilities:
    u: np.ndarray
    v: np.ndarray
    weight : float 
    V_reals : np.ndarray
    V_imags : np.ndarray
    nchan : int 
    nbaselines : int

    def __post_init__(self):
        self.V = self.V_reals + 1.0j * self.V_imags

        # initialize Hermitian pairs
        self.uu = np.concatenate([self.u, -self.u]) # shape (2 * nsamples)
        self.vv = np.concatenate([self.v, -self.v]) # shape (2 * nsamples)
        self.VV = np.concatenate([self.V, np.conj(self.V)], axis=1) # shape (nchan, 2 * nsamples)
        self.weightweight = self.weight * np.ones_like(self.uu) # shape (2 * nsamples)
        self.weight = self.weight * np.ones_like(self.u)

nu0 = 345.79598990e9 # [Hz]
lam0 = c.value / nu0  # m

# load the uncorrelated dataset
d = np.load("uncorrelated_dataset.npz")

u = d["u"]
v = d["v"]

# convert to klambda
u = 1e-3 * u / lam0  # [klambda]
v = 1e-3 * v / lam0  # [klambda]

# assume channel spacings are small enough such that baseline values are the same in adjacent channels

nchan, nbaselines = d["V_imags"].shape

vis_uncorrelated = Visibilities(u, v, d["weight"], d["V_reals"], d["V_imags"], nchan, nbaselines)



# load the correlated dataset
d = np.load("correlated_dataset.npz")

u = d["u"]
v = d["v"]

# convert to klambda
u = 1e-3 * u / lam0  # [klambda]
v = 1e-3 * v / lam0  # [klambda]

# assume channel spacings are small enough such that baseline values are the same in adjacent channels

nchan, nbaselines = d["V_imags"].shape

vis_correlated = Visibilities(u, v, d["weight"], d["V_reals"], d["V_imags"], nchan, nbaselines)
