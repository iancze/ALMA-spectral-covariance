import numpy as np

if __name__ == "__main__":
    import matplotlib.pyplot as plt 
    import casatools

    msmd = casatools.msmetadata()
    spw_id = 0
    msmd.open("LkCa15-spw3.ms")
    # TOPO
    chan_freq = msmd.chanfreqs(spw_id)
    msmd.done()

    nchan = len(chan_freq)

    window = np.hanning(nchan + 1)[:-1]

    N_factor_pad = 10
    window_packed = np.concatenate((np.zeros(N_factor_pad//2 * len(window)), window, np.zeros(N_factor_pad//2 * len(window))))
    nchan_padded = len(window_packed)


    response = np.fft.fftshift(np.fft.fft(np.fft.fftshift(window_packed)))



    fig, ax = plt.subplots(nrows=3)
    
    ax[0].plot(window_packed)
    
    ax[1].plot(response.real)
    ax[1].plot(response.real, ".")
    ax[2].plot(response.imag)

    for i in [1,2]:
        ax[i].set_xlim(nchan_padded//2 - 5 * N_factor_pad, nchan_padded//2 + 5 * N_factor_pad)
    plt.savefig("apodization-functions.png", dpi=200)