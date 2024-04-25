import numpy as np

# https://github.com/astroML/astroML/blob/main/astroML/fourier.py

def sinegauss(t, t0, f0, Q):
    """Sine-gaussian wavelet. Complex-valued."""
    a = (f0 * 1. / Q) ** 2
    return (np.exp(-a * (t - t0) ** 2)
            * np.exp(2j * np.pi * f0 * (t - t0)))

def cosgauss(t, t0, f0, sigma):
    """A Gaussian modified by a cosine window. Entirely real.
    
    Args:
        t: time coordinate
        t0: central time coordinate, where Gaussian is centered
        f0: frequency of cosine modulation.
        sigma: width of Gaussian.
        
    Returns:
        cosine Gauss evaluated at t"""

    return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-(t - t0)**2 / (2 * sigma**2)) * np.cos(2 * np.pi * f0 * (t - t0))
    
def raisedcosgauss(x, x0, f0, sigma):
    """Similar to cosgauss but with an extra added Gaussian bump. That way,
    the function is entirely real and non-negative. The minimum value is 0.
    
    Args:
        t: central x coordinate
        x0: central x coordinate, where Gaussian is centered
        f0: frequency of cosine modulation.
        sigma: width of Gaussian.
        
    Returns:
        cosine Gauss evaluated at x"""

    return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - x0)**2 / (2 * sigma**2)) * (np.cos(2 * np.pi * f0 * (x - x0)) + 1)
    

def sinegauss_FT(f, t0, f0, Q):
    """Fourier transform of the sine-gaussian wavelet.
    This uses the convention
    .. math::
       H(f) = integral[ h(t) exp(-2pi i f t) dt]
    """
    a = (f0 * 1. / Q) ** 2
    return (np.sqrt(np.pi / a)
            * np.exp(-2j * np.pi * f * t0)
            * np.exp(-np.pi ** 2 * (f - f0) ** 2 / a))

def cosgauss_FT(f, t0, f0, sigma):
    """Fourier transform of the sine-gaussian wavelet.
    This uses the convention
    .. math::
       H(f) = integral[ h(t) exp(-2pi i f t) dt]
    """
    return np.exp(-2.0j * np.pi * t0 * f)/2 * (np.exp(-2 * np.pi**2 * sigma**2 * (f - f0)**2)  + np.exp(-2 * np.pi**2 * sigma**2 * (f + f0)**2))

def raisedcosgauss_FT(f, t0, f0, sigma):
    """Fourier transform of the raised Gauss - cosine window.
    """
    return np.exp(-2.0j * np.pi * t0 * f)/2 * (np.exp(-2 * np.pi**2 * sigma**2 * (f - f0)**2)  + np.exp(-2 * np.pi**2 * sigma**2 * (f + f0)**2) + 2 * np.exp(-2 * np.pi**2 * sigma**2 * f**2))

def fftshift(xs, vals, nu_shift):
    """
    Shift the function by applying a phase shift in the frequency domain

    Args:
        nu_shift: how much to shift by, in units of Hz
    """
    N = len(vals)
    assert N % 2 == 0, "Only even numbered arrays for now."
    dchan = xs[1] - xs[0]

    # use fft to access cross-correlation function
    rho_packed = np.fft.fft(np.fft.fftshift(vals))
    fs_packed = np.fft.fftfreq(N, d=dchan)

    # mulitply rho_packed by phase shift 
    # if shifting by positive dnu, then 
    phase = np.exp(-2.0j * np.pi * fs_packed * nu_shift)

    # transform back using ifft
    rho_packed_shifted = rho_packed * phase 
    return np.fft.fftshift(np.fft.ifft(rho_packed_shifted))

def fftshift_data(xs, vals, nu_shift):
    """
    Shift the function by applying a phase shift in the frequency domain

    Args:
        nu_shift: how much to shift by, in units of Hz
    """
    N = vals.shape[1]
    assert N % 2 == 0, "Only even numbered arrays for now."
    dchan = xs[1] - xs[0]

    # use fft to access cross-correlation function
    rho_packed = np.fft.fft(np.fft.fftshift(vals, axes=1), axis=1)
    fs_packed = np.fft.fftfreq(N, d=dchan)

    # mulitply rho_packed by phase shift 
    # if shifting by positive dnu, then 
    phase = np.exp(-2.0j * np.pi * fs_packed * nu_shift)

    # transform back using ifft
    rho_packed_shifted = rho_packed * phase[np.newaxis, :, np.newaxis] 
    return np.fft.fftshift(np.fft.ifft(rho_packed_shifted, axis=1), axes=1)




def apodize(vals):
    """
    Take a wavelet, then FFT, apply Hann window, then take back.
    """
    nchan = len(vals)
    rho_packed = np.fft.fft(np.fft.fftshift(vals))
    window = np.hanning(nchan + 1)[:-1]
    window_packed = np.fft.fftshift(window)
    rho_packed_windowed = window_packed * rho_packed 

    return np.fft.fftshift(np.fft.ifft(rho_packed_windowed))

    
def apodize_fine(chan_freq, vals, N_factor_pad=10):
    """
    Take a wavelet, FFT, apply Hann window, then back.

    But, also pad with zeros by a large factor so we can see the SRF.
    """

    nchan = len(chan_freq)
    window = np.hanning(nchan + 1)[:-1]
    
    rho = np.fft.fftshift(np.fft.fft(np.fft.fftshift(vals)))
    rho_windowed = window * rho
    
    rho_window_zeros = np.concatenate((np.zeros(N_factor_pad//2 * len(window)), rho_windowed, np.zeros(N_factor_pad//2 * len(window))))

    # multiply by a pre-factor to keep normalization the same 
    prefactor = N_factor_pad + 1
    rho_window_zeros_packed = prefactor * np.fft.fftshift(rho_window_zeros)
    
    vals_windowed_sampled = np.fft.fftshift(np.fft.ifft(rho_window_zeros_packed))

    # calculate the new channel frequencies 
    # there should now be N_factor_pad additional channel points between each channel

    dchan_orig = chan_freq[1] - chan_freq[0]
    dchan_new = dchan_orig / (N_factor_pad + 1)

    nchan_new = (N_factor_pad + 1) * nchan

    chan_freq_new = chan_freq[0] + np.arange(nchan_new) * dchan_new
    
    return chan_freq_new, vals_windowed_sampled



if __name__ == "__main__":
    import matplotlib.pyplot as plt 
    import casatools
    import barytimes
    import visread.utils

    msmd = casatools.msmetadata()

    spw_id = 0
    msmd.open("LkCa15-spw3.ms")

    # TOPO
    chan_freq = msmd.chanfreqs(spw_id)
    msmd.done()

    # fftshift works by shifting the phase of the signal, so a linear
    # shift in channel number. What is the delta v implied by delta 
    # nu at the start of the bandpass, and delta v at the end of the 
    # bandpass? 
    # vstart = (chan_freq[1] - chan_freq[0]) / chan_freq[0]
    # vend = (chan_freq[-1] - chan_freq[-2]) / chan_freq[-1]

    # print(vstart, vend)
    # print((vstart - vend)/vstart)
    # # 0.01% of channel width, so we should be fine.


    # assume that chan_freq is TOPO frame.
    # what channel frequencies would these correspond to in BARY frame
    # convention is that barytimes is *added* to velocity
    # so, if we assume that chan_freq corresponds to BARY
    chan_freq_epoch = visread.utils.doppler_shift(chan_freq, barytimes.corrs[0])

    print("delta nu bary", chan_freq_epoch[0] - chan_freq[0])

    # print("bary 0", barytimes.corrs[0])
    # - 10.26 km/s to go from TOPO to BARY

    N = 3840
    dchan = chan_freq_epoch[1] - chan_freq_epoch[0]
    
    # BARY frame coordinates
    x0 = 345.765e9 # Hz
    f0 = 1.7e-05 # oscillations per Hz
    sigma = 30000 # Hz

    # should be centered time-stream on 0
    chan_number = (chan_freq_epoch - np.min(chan_freq_epoch)) / dchan

    xs_fine = np.linspace(np.min(chan_freq_epoch), np.max(chan_freq_epoch), num=100000)
    hs_fine = raisedcosgauss(xs_fine, x0, f0, sigma)

    print("central channel", (x0 - np.min(chan_freq_epoch))/dchan) # 2021
    
    # set two x-axes to represent channel number and frequency

    hs = raisedcosgauss(chan_freq_epoch, x0, f0, sigma)
    hs_hann = apodize(hs)
    xs_hann_fine, hs_hann_fine = apodize_fine(chan_freq_epoch, hs)

    fig, ax = plt.subplots(nrows=3, sharex=True)
    ax[0].plot(xs_fine, hs_fine.real)
    ax[0].plot(xs_hann_fine, hs_hann_fine.real, label="Hann")
    ax[0].plot(chan_freq_epoch, hs.real, ".", "rect")
    ax[0].plot(chan_freq_epoch, hs_hann.real, ".", label="Hann")
    ax2 = ax[0].twiny()

    ax[0].legend()

    shifted_hs = fftshift(chan_freq_epoch, hs, 1.5 * dchan)
    hs_fine_shifted = raisedcosgauss(xs_fine, x0 + 1.5 * dchan, f0, sigma)
    ax[0].plot(xs_fine, hs_fine_shifted.real)
    ax[0].plot(chan_freq_epoch, shifted_hs.real, ".")   

    ax[1].plot(xs_fine, hs_fine.imag)
    ax[1].plot(chan_freq_epoch, hs.imag, ".")
    ax[1].plot(chan_freq_epoch, hs_hann.imag, ".", label="Hann")
    ax[1].plot(chan_freq_epoch, shifted_hs.imag, ".")   
    ax[1].legend()

    ax[2].plot(xs_fine, np.abs(hs_fine))
    ax[2].plot(chan_freq_epoch, np.abs(hs), ".")
    ax[2].plot(chan_freq_epoch, np.abs(hs_hann), ".")
    ax[2].set_xlim(x0 - 10 * sigma, x0 + 10 * sigma)
    ax2.set_xlim(((x0 - 10 * sigma) - np.min(chan_freq_epoch))/dchan, ((x0 + 10 * sigma) - np.min(chan_freq_epoch))/dchan)
    ax2.set_xlabel("channel number")
    
    ax[0].set_ylabel(r"$\Re f$")
    ax[1].set_ylabel(r"$\Im f$")
    ax[2].set_ylabel(r"$|f|$")
    ax[2].set_xlabel(r"$\nu$ [BARY, Hz]")
    fig.savefig("wavelet-function.png", dpi=200)

    # now compute the FT and compare the analytic version
    Hs_FFT = dchan * np.fft.fftshift(np.fft.fft(np.fft.fftshift(hs)))
    freqs = np.fft.fftshift(np.fft.fftfreq(N, d=dchan))
    Hs_analytical = raisedcosgauss_FT(freqs, x0, f0, sigma)

    fig, ax = plt.subplots(nrows=2, sharex=True)
    ax[0].plot(freqs, Hs_analytical.real)
    ax[0].plot(freqs, Hs_FFT.real, ".")
    
    ax[1].plot(freqs, Hs_analytical.imag)
    ax[1].plot(freqs, Hs_FFT.imag, ".")

    ax[0].set_ylabel(r"$\Re F$")
    ax[1].set_ylabel(r"$\Im F$")
    ax[1].set_xlabel(r"$s$")
    fig.savefig("wavelet-fft-function.png", dpi=200)
