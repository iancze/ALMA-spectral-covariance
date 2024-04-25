# make a gif stepping through the channels
import numpy as np 
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt 
from matplotlib.animation import FFMpegWriter
from tqdm import tqdm

metadata = dict(title='Uncorrelated channels', artist='Ian Czekala',
                comment='Each synthesized image contains only noise.')

writer = FFMpegWriter(fps=1, metadata=metadata, bitrate=20000, codec="mpeg4")


# load the uncorrelated cube and set limits from that. Re-use limits for Hann, since they should be smaller.
dcube = np.load("image_cube_uncorrelated.npz")
img_cube = dcube["img_cube"]
# get the min/max intensities across the cube
vmin = np.min(img_cube)
vmax = np.max(img_cube)

def plot_image_gif(cube_name, plot_name):
    
    dcube = np.load(cube_name)
    img_cube = dcube["img_cube"]
    ls = dcube["ls"]
    ms = dcube["ms"]
    nchan = len(img_cube)

    
    dl = np.abs(ls[1] - ls[0])

    # (left, right, bottom, top),
    extent = [np.max(ls) + 0.5 * dl, np.min(ls) - 0.5 * dl, np.min(ms) - 0.5 * dl, np.max(ms) + 0.5 * dl]

    xx = 4 # in
    lmargin = 0.5 # in
    rmargin = 0.1 # in
    bmargin = 0.5 # in
    tmargin = 0.02 # in 

    ax_width = (xx - lmargin - rmargin) # in
    ax_height = ax_width # in

    yy = bmargin + ax_height + tmargin

    fig = plt.figure(figsize=(xx, yy))
    
    # (left, bottom, width, height)
    a = fig.add_axes([lmargin/ xx, bmargin/yy, ax_width/xx, ax_height/yy])
    a.tick_params(direction="in")

    a.set_xlabel(r"$\Delta \alpha \cos \delta$ [arcseconds]")
    a.set_ylabel(r"$\Delta \delta$ [arcseconds]")

    
    im = a.imshow(np.fliplr(img_cube[0]), origin="lower", extent=extent, vmin=vmin, vmax=vmax, cmap="inferno")
    # channel label 
    t = a.annotate("chan = {:2d}".format(0), (-1.3, 2.7), color="w", fontsize="large")
     
    with writer.saving(fig, plot_name, 300):
        # writer.grab_frame()
        for i in tqdm(range(nchan)):
            im.set_data(np.fliplr(img_cube[i]))
            t.set_text("chan = {:2d}".format(i))
            writer.grab_frame()


plot_image_gif("image_cube_uncorrelated.npz", "image_uncorrelated.mp4")
plot_image_gif("image_cube_correlated.npz", "image_correlated.mp4")