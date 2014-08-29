RRLspectra
==========

IDL libraries for processing low frequency CRRL spectra

To process spectra, simply run process_spectra.pro

NB:
1) This code makes use of the coyote IDL libraries (http://idlcoyote.com/documents/programs.php)<\br>
2) A hard-coded path to an LLINE directory (included in this download) is required
3) start from a new directory that has all of the *.restored.corr images
4) there are two options for making the fits image: with smoothing, or without smoothing
   NO SMOOTHING: makefits.py
   SMOOTHING: makefits_sm.py
   --> you will have to change this (line X in process_spectra.pro) and give the path to the python script
   --> The smoothing script smooths to 400x400 arcsec. You can change this value in the script if you want.
   --> If you run without smoothing first, process_spectra will tell you what the smallest and largest beam sizes are
   --> and then you can run again with a kernal that is slightly larger than the largest beam size!

process_spectra will then do everything else for you, except for stacking (which requires an aperture size and a redshift)

To stack, run stack.pro:
stack, sigma, redshift
This is meant to be run from the same directory that you ran process_spectra in, and it will look for specific directories
created by process_spectra

Scripts that process_spectra utilizes:
  - makefits.py (or makefits_sm.py)
  - aperextraction.pro
  - apersbplot.pro

These will need to be in your IDL path (ideally) or in the directory you run process_spectra from.

For questions contact L.Morabito at morabito@strw.leidenuniv.nl
