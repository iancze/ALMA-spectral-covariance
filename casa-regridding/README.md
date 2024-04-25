# Spectral correlations

11/28/22

We investigated how spectral correlations worked in a simple example.

Now, let's use CASA and a fake measurement set, on the scale of the exoALMA data, to demonstrate the problem exactly.

This may lead us to implementing an fftshift interpolation capability, or not.


- [x] Check to ensure wavelet and FFT are implemented consistently
- [x] Split out short-baseline LkCa 15 dataset in TOPO frame 
- [x] Image a few channels to get an idea of the noise
- [x] Replace visibilities with wavelet function centered on channels, in TOPO frame
- [x] Image with cubedata to see what tclean produces, hopefully the amplitude of the point source matches the expected spectral profile. This didn't exactly work, because I think we need to be more careful about how we are inserting real and imaginary components. We have assigned the point source a complex spectrum, which is OK *if* the point source were not at the phase center, or it were an extended source. However, if we want this point source to be at center, then we need to assign it a truly real spectrum.

We found that CASA > 6.2 has a bug in cubedata, so we can't use this. However, we verified that our routine based on a single integration works fine.

- [x] To correctly insert the profile, we need to calculate what the barycentric correction would have been for each measurement, and then insert that into the data as if it were Doppler tracking.
- [x] Image the raw dataset using tclean in the LSRK frame, using linear interpolation, compare spectral profile. We found that the profile was substantially blurred. Finer features removed.
- [x] Examine spectral profile using cubic regridding. Not substantially changed from linear.

Overall, though, I have been experiencing many issues with the tclean approach. The main thing is that there is no 'cubedata' option that works in casa > 6.2, so I can't directly verify that the data is as I would expect. I think it would be better if I just did everything in cvel2, and then I could show the 1D profiles. I'm getting some strange behavior when I image in TOPO, with no-regridding. I think we verified that tclean messes things up when linear interpolation is used to regrid channels. 

From the cvel2 approach, what would we want?

- [x] refactor wavelet injection routine to focus on channel frequencies
- [x] insert wavelet into an MS via Doppler setting

* there are 65 individual times in this measurement set. Can we show the individual interpolations? Problem is we'd need to show them for all 741 baselines, too. But, I think the point is to show the samples are different for each slice.
- [x] use cvel2 and linear regridding to produce an MS with all visibilities in BARY frame. Also try cubic.

* show average of linearly interpolated values, as well as scatter
* show average of fftshift interpolated values


Do I need to worry about the spectral response function of the correlator? I just put the full-resolution cosine-Gaussian wavelet into the visibility domain. I think this is OK because the function was already band-limited. 

But, we will need to consider how the Hann window affects this regridding operation to make a comparison.

So, technically, I should add the FT of the cos-Gauss window to the cross-correlation function, modified by a phase shift, and then modified by a Hann window. 
Then, FFT back and fill up the frequency spectrum. That way, this will have the SRF imprints.

* let's make a plot comparing the Hann wavelet to the raw wavelet

And, it will be easy enough to implement the FFT shift back.

I think the answer is that this shouldn't change for the sinc/boxcar limit, because the signal is already band-limited. The FFT and agreement of the samples shows this, but it would be even better if we could show it numerically (fft - true).

* plot regridded linear spectrum relative to what we know the injected spectrum to be
* demonstrate that we could have recovered the true injected spectrum with fftshift capability.
* is the difference between the raw spectrum and that produced via linear interpolation a convolution with a triangle kernel?

- [ ] Bin to proposed exoALMA channel widths and examine spectral profile. Compare to 
- [ ] Implement fftshift capability and see if result still matches TOPO frame result.


So we proved that we can exactly recover the signal. But, we are worried that there are actual errors in the different barycentric routine corrections. It's clear the first transform has a 5 m/s error. How does this change over the duration of the EB?
