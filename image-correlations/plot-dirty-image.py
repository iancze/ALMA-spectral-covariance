import numpy as np 
import matplotlib.pyplot as plt 

# load the uncorrelated cube and set limits from that. Re-use limits for Hann, since they should be smaller.
dcube = np.load("image_cube_uncorrelated.npz")
img_cube = dcube["img_cube"]
# get the min/max intensities across the cube
# get the min/max intensities across the channels of interest 

# Choose a channel to image
nchan = len(img_cube)
chan0 = nchan//2 - 3

vmin = np.min(img_cube[chan0:chan0+5])
vmax = np.max(img_cube[chan0:chan0+5])

def plot_image_cube(cube_name, plot_name):
    
    dcube = np.load(cube_name)
    img_cube = dcube["img_cube"]
    ls = dcube["ls"]
    ms = dcube["ms"]
    nchan = len(img_cube)


    # create a 5-panel figure
    npanels = 5

    # vmin = np.min(img_cube[chan0:chan0+5])
    # vmax = np.max(img_cube[chan0:chan0+5])

    dl = np.abs(ls[1] - ls[0])

    # (left, right, bottom, top),
    extent = [np.max(ls) + 0.5 * dl, np.min(ls) - 0.5 * dl, np.min(ms) - 0.5 * dl, np.max(ms) + 0.5 * dl]

    xx = 11 # in
    lmargin = 0.5 # in
    rmargin = 0.5 # in
    bmargin = 0.5 # in
    tmargin = 0.3 # in 
    mmargin = 0.05 # in

    ax_width = (xx - lmargin - rmargin - (npanels - 1) * mmargin) / npanels # in
    ax_height = ax_width # in

    yy = bmargin + ax_height + tmargin

    fig = plt.figure(figsize=(xx, yy))
    ax = []
    for i in range(npanels):
        chan_ind = chan0 + i
        # (left, bottom, width, height)
        a = fig.add_axes([(lmargin + (ax_width + mmargin) * i) / xx, bmargin/yy, ax_width/xx, ax_height/yy])
        im = a.imshow(np.fliplr(img_cube[i]), origin="lower", extent=extent, vmin=vmin, vmax=vmax, cmap="inferno")
        a.tick_params(direction="in")
        ax.append(a)

    for a in ax[1:]:
        a.xaxis.set_ticklabels([])
        a.yaxis.set_ticklabels([])
        # a.axis("off")

    # set titles
    ax[0].set_title(r"$n - 2$")
    ax[1].set_title(r"$n - 1$")
    ax[2].set_title(r"$n$")
    ax[3].set_title(r"$n + 1$")
    ax[4].set_title(r"$n + 2$")

    # ax_cb = plt.colorbar(im, ax=ax[-1])
    # ax_cb.set_label(r"$I$ [Jy/dirty beam]")
    ax[0].set_xlabel(r"$\Delta \alpha \cos \delta$ [arcseconds]")
    ax[0].set_ylabel(r"$\Delta \delta$ [arcseconds]")

    fig.savefig(plot_name, dpi=200)

plot_image_cube("image_cube_uncorrelated.npz", "image_uncorrelated.png")
plot_image_cube("image_cube_correlated.npz", "image_correlated.png")