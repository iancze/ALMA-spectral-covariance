import numpy as np 
import matplotlib.pyplot as plt 



boxcar_width = 4.0

# stretches the full boxcar
N = 40
dtau = 0.1 
taus = dtau * np.arange(-N/2, N/2)

# goes further than boxcar, and is more finely sampled
N_fine = N * 20
dtau_fine = dtau * (N/(N_fine - 8 * N))

taus_fine = dtau_fine * np.arange(-N_fine/2, N_fine/2)

nu_0 = 1.5

def sine_func(tau):
    return np.sin(2* np.pi * nu_0 * tau)



@np.vectorize
def boxcar(tau):
    if np.abs(tau) < boxcar_width/2:
        return 1.0
    else:
        return 0.0

rhos = sine_func(taus)
rhos_fine = sine_func(taus_fine)

boxcar_fine = boxcar(taus_fine)

window = np.hanning(N + 1)[:-1]
window_fine = np.hanning(N_fine + 1)[:-1]

def hann_extended(tau):
    # actual functional form
    hann = 1/2 + 0.5 * np.cos(2 * np.pi * tau / boxcar_width)
    return boxcar(tau) * hann
    

fig, ax = plt.subplots(nrows=3, figsize=(6, 6))
ax[0].plot(taus_fine, rhos_fine)
# ax[0].plot(taus, rhos, ".")

ax[1].plot(taus_fine, rhos_fine * boxcar_fine)
ax[1].plot(taus_fine, boxcar_fine)

ax[2].plot(taus_fine, rhos_fine * hann_extended(taus_fine))
ax[2].plot(taus_fine, hann_extended(taus_fine))

for a in ax:
    a.axis("off")

fig.savefig("taus.png", dpi=300)

freqs = np.fft.fftshift(np.fft.fftfreq(N, d=dtau))

rhos_rect_windowed = boxcar(taus) * rhos
spectrum = dtau * np.fft.fftshift(np.fft.fft(np.fft.fftshift(rhos_rect_windowed)))


# pack the rhos with zeros to get more points
N_factor_pad = 10
rhos_rect_windowed_padded = np.concatenate((np.zeros(N_factor_pad//2 * len(rhos_rect_windowed)), rhos_rect_windowed, np.zeros(N_factor_pad//2 * len(rhos_rect_windowed))))

window = np.hanning(N + 1)[:-1]
rhos_hann_windowed = window * rhos
rhos_hann_window_padded = np.concatenate((np.zeros(N_factor_pad//2 * len(rhos_hann_windowed)), rhos_hann_windowed, np.zeros(N_factor_pad//2 * len(rhos_hann_windowed))))
spectrum_hann = dtau * np.fft.fftshift(np.fft.fft(np.fft.fftshift(rhos_hann_windowed)))

dchan_orig = freqs[1] - freqs[0]
dchan_new = dchan_orig / (N_factor_pad + 1)
nchan_new = (N_factor_pad + 1) * N
freqs_fine = freqs[0] + np.arange(nchan_new) * dchan_new
    
spectrum_rect_fine = dtau * np.fft.fftshift(np.fft.fft(np.fft.fftshift(rhos_rect_windowed_padded)))
spectrum_hann_fine = dtau * np.fft.fftshift(np.fft.fft(np.fft.fftshift(rhos_hann_window_padded)))

fig, ax = plt.subplots(nrows=3, figsize=(6, 6), sharex=True)

arrow_props = {"width": 0.01, "head_width": 0.1, "color":"k"}
ax[0].axhline(0, c='k', lw=1)
ax[0].arrow(-nu_0, 0, 0, 2, **arrow_props)
ax[0].arrow(nu_0, 0, 0, -2, **arrow_props)
ax[0].set_ylim(-2.5, 2.5)

ax[1].plot(freqs_fine, spectrum_rect_fine.imag, c="C0", label="rect")
ax[1].plot(freqs_fine, spectrum_hann_fine.imag, c="C1", label="Hann")
ax[1].plot(freqs, spectrum.imag, ".", c="C0")
ax[1].plot(freqs, spectrum_hann.imag, ".", c="C1")
ax[1].legend()

ax[2].plot(freqs_fine, np.abs(spectrum_rect_fine), c="C0")
ax[2].plot(freqs_fine, np.abs(spectrum_hann_fine), c="C1")
ax[2].plot(freqs, np.abs(spectrum), ".", c="C0")
ax[2].plot(freqs, np.abs(spectrum_hann.imag), ".", c="C1")

ax[2].set_xlim(np.min(freqs), np.max(freqs))

fig.savefig("nus.png", dpi=300)


# def plot_fig(xs, ys, fname):
#     fig, ax = plt.figure()
