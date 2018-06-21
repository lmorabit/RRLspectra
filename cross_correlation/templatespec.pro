function templatespec, fmin, fmax, c_freq, vwidth

	;; constants
	speedoflight = 2.99792458d5 ;; km/s

;	if trans eq 'alpha' then readcol,'LLINE/RRL_CI_alpha.txt',ca1, ca2, cann, c_freq, format='A,A,D,D'
;	if trans eq 'beta' then readcol,'LLINE/RRL_CI_beta.txt',cb1, cb2, cbnn, c_freq, format='A,A,D,D'
;	if trans eq 'gamma' then readcol,'LLINE/RRL_CI_gamma.txt',cg1, cg2, cgnn, c_freq, format='A,A,D,D'

	cfreqs = c_freq[where(c_freq ge fmin and c_freq le fmax)] ;; already in MHz
	cfreqs = reverse(cfreqs) ;; to use value locate, numbers must be monotonically increasing

	;; convert to frequency
	gwidth = cfreqs * vwidth / speedoflight
	;; and translate to gaussian sigma
	gsig = gwidth / 2d0 / sqrt(2d0*alog(2d0))

	;; define a continuum level
	nabs = (fmax - fmin)*1d5
	cont_abscissa = dindgen(nabs)/1d5 + fmin
	cont_values = dblarr(n_elements(cont_abscissa))

	;; find where cfreqs are in the abscissa
	cfreqs_index = Value_Locate(cont_abscissa, cfreqs)

	;; set up gaussian parameters
	nvals = 10001
	xvals = dindgen(nvals)/1d5 - double(nvals/2)/1d5
	maxval = 1d0
	
	for ii=0,n_elements(cfreqs_index)-1 do begin

		gs = gsig[ii]
		tmpg = gaussian(xvals, [maxval, 0d0, gs])
		cont_values[cfreqs_index[ii]-(nvals/2):cfreqs_index[ii]+(nvals/2)] = cont_values[cfreqs_index[ii]-(nvals/2):cfreqs_index[ii]+(nvals/2)] - tmpg
	
	endfor ;; index ii


	tempspec = dblarr(2,n_elements(cont_abscissa))
	tempspec[0,*] = cont_abscissa
	tempspec[1,*] = cont_values

return, tempspec
end
	
