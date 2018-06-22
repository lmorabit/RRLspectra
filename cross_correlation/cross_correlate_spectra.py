import numpy as np
import argparse
import pandas
import matplotlib.pyplot as plt
import time
from scipy.interpolate import interp1d

def find_nearest(array, value):
	idx = (np.abs(array - value)).argmin()
	return idx

def find_array(myArray, refArray):
	out1 = np.empty(myArray.size, dtype=int)
	for i, value in enumerate(myArray):
        	index = find_nearest(refArray, value)
		out1[i] = index
	return out1

def templatespec(fmin, fmax, c_freq, vwidth):
	speedoflight = 2.99792458e5

	cfreqs = [i for i in c_freq if i >= fmin and i <=fmax]
	cfreqs = np.asarray(cfreqs[::-1])
	gwidth = cfreqs * vwidth / speedoflight
	gsig = gwidth / 2.0 / np.sqrt(2.0* np.log(2.0))

	nabs = (fmax -fmin) * 1e5
	cont_abscissa = np.arange(nabs)/1e5 + fmin
	cont_values = np.zeros(cont_abscissa.size)

	#print('cont_abscissa',cont_abscissa.size)
	cfreqs_index = find_array(cfreqs,cont_abscissa)

	#print('cfreqs_index',cfreqs_index)

	nvals = 10000
	xvals = np.arange(nvals) / 1e5 - nvals/ 2 / 1e5
	maxval = 1e0

	for ii in range(0,len(cfreqs_index)-1):

		gs = gsig[ii]
		tmpg = maxval*np.exp(-(xvals-0.0)**2/gs**2)
	        cont_values[cfreqs_index[ii] - int((nvals / 2)):cfreqs_index[ii] + int((nvals / 2))] = cont_values[cfreqs_index[ii] - int((nvals / 2)):cfreqs_index[ii] + int((nvals / 2))] - tmpg

	tempspec = np.zeros((2, cont_abscissa.size))
	tempspec[0, :] = cont_abscissa
	tempspec[1, :] = cont_values

	return tempspec

def get_observed_spectrum( filename ):

	obs_spec_array = np.loadtxt( filename )
	obs_freq = obs_spec_array[:,0]
	obs_opt_depth = obs_spec_array[:,1] / 100.  ## convert from percent to values

	return obs_freq, obs_opt_depth


def main( spectral_file, vel_min=0., vel_max=0., vel_step=0., zmin=0., zmax=0., zstep=0., crrl_file='/home/morabito/scripts/RRLspectra/cross_correlation/LLINE/RRL_CI_alpha.txt' ):

        ## get observed spectrum
        obs_freq, obs_opt_depth = get_observed_spectrum( spectral_file )
	fmin = np.min( obs_freq )
	fmax = np.max( obs_freq )

	## read in rest frequencies of lines
	crrl_table = np.genfromtxt(crrl_file)
	crrl_freqs = crrl_table[:,3]

	## define the velocity widths / step
	vel_widths = np.arange( vel_min, vel_max, vel_step )
	## define the redshift range /step
	redshifts = np.arange( zmin, zmax, zstep )

	## set up an array of len( vel_widths ) x len( redshifts )
	cross_corr_coeffs = np.zeros((len(vel_widths),len(redshifts)))
	
	## loop over velocity width
	for ii in np.arange(0,len(vel_widths)):
		print( 'vel width', vel_widths[ii] )
		## get the template spectrum for the velocity width
		template_spectrum = templatespec( fmin, fmax*(1.+max(redshifts)), crrl_freqs, vel_widths[ii] ) #TA, redshifting the template takes it out of range
		template_freq = template_spectrum[0,:]
		template_opt_depth = template_spectrum[1,:]

		## loop over redshift
		for jj in np.arange(0,len(redshifts)):
			print( '... redshift', redshifts[jj] )
			## redshift the template spectrum and sample the same as the observed
			zz = redshifts[jj]
			template_freq_obs = template_freq / ( 1. + zz )

			## find indices of where there is observed coverage
			interp_func = interp1d( template_freq_obs, template_opt_depth, kind='nearest', bounds_error=False, fill_value='extrapolate' )
			sampled_template_opt_depth = interp_func( obs_freq )

			#print('indexing...')
			#t1=time.time()
			#index = find_array( obs_freq, template_freq_obs )
			#sampled_template_opt_depth = template_opt_depth[index]
			#t2=time.time()
			#print('finished indexing, time='t2-t1,' s')

			## drop all the nans
			finite_index = np.where( np.isfinite( obs_opt_depth ) )

			observed_opt_depth = obs_opt_depth[finite_index]
			samp_template_opt_depth = sampled_template_opt_depth[finite_index]
	
			## cross-correlate
			print('running cross-correlation...')
			t3=time.time()
			print( samp_template_opt_depth[1:20] )
			print( observed_opt_depth[1:20] )
			cross_corr_coeffs[ii,jj] = np.corrcoef( samp_template_opt_depth, observed_opt_depth )[1,0]
			t4=time.time()
			print('finished cross-correlation! time=',t4-t3,' s')

	outfile = spectral_file.replace('.txt','_xcorr.txt')
	df = pandas.DataFrame(cross_corr_coeffs, columns=redshifts, index=vel_widths)
	df.to_csv( outfile )


if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument( '--vmin', dest='vel_min', type=float, help='minimum velocity width [km/s] of line' )
	parser.add_argument( '--vmax', dest='vel_max', type=float, help='minimum velocity width [km/s] of line' )
	parser.add_argument( '--vstep', dest='vel_step', type=float, default=2., help='step for velocity width [km/s]' )
	parser.add_argument( '--zmin', dest='zmin', type=float, help='minimum redshift' )
	parser.add_argument( '--zmax', dest='zmax', type=float, help='maximum redshift' )
	parser.add_argument( '--zstep', dest='zstep', type=float, default=0.1, help='step for redshift' )
	parser.add_argument( 'spectral_file', type=str, help='filename for spectrum to search' )
	
	args = parser.parse_args()
	main( args.spectral_file, vel_min=args.vel_min, vel_max=args.vel_max, vel_step=args.vel_step, zmin=args.zmin, zmax=args.zmax, zstep=args.zstep )
