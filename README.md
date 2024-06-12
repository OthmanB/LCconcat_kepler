# LCconcat_kepler
 Create concatenated lightcurves from Kepler MAST data

Updated on 12 June 2024: 
The Python functions to create the LombScargle periodogram were either slow or not woorking as I expected.
I therefore reversed back to a computation of the PSF using the IDL code. It is far from ideal, but it works
perfectly. 
If someday I figure out how to get the same outputs as the IDL lnp_test() function, I will get it implemented.
