pro ncrosscorr2
;
;
;	Cross correlate a template spectrum of RRLs (alpha, beta, or gamma transitions) with an observation
;
;	Reads in template and observed spectra, cross correlates, outputs ....
;
;	Written by Leah Morabito
;	version 1.0 01/23/2013: Initial version

	set_plot,'ps'
	!P.font=0

	;; read in observed spectrum
	sigma = 2.2
	sigmastr = cgnumber_formatter(sigma,decimal=1)
	specfile = '/net/leubeek/data1/m82/sgstack/TOTALSPEC_' + sigmastr + '.dat'
	readcol, specfile, frequency, spectra,weights, subband, format='D,D,D,D'

	;; read in frequencies for template spectrum
	readcol,'/net/leubeek/data1/m82/LLINE/RRL_CI_alpha.txt',ca1, ca2, cann, c_freq, format='A,A,D,D'

	;; clear the nans
	freqs = frequency[where(finite(frequency) ne 0)] 
	spec = spectra[where(finite(spectra) ne 0)]
	chwidth = freqs[1] - freqs[0]

	;; construct fake subband identifiers
	nsbs = n_elements(frequency) - n_elements(freqs)
	nanindex = where(finite(frequency) eq 0)
	nchans = nanindex[1] - nanindex[0]-1
	sbs = indgen(nsbs)+1
	subbands = replicate(1,nchans)
	for ii=1,n_elements(sbs)-1 do subbands = [subbands, replicate(sbs[ii],nchans)]

	freqmin = min(freqs)
	freqmax = max(freqs)

	;; redshift range
        ;redshiftrange = (dindgen(31)+60d0)/1d5 ;; redshift from 0.0000 to 0.0009 in steps of 0.00001
        redshiftrange = (dindgen(101))/1d5 ;; redshift from 0.0000 to 0.0009 in steps of 0.00001

	;; different line widths in km/s
	linevelocities = dindgen(25)+5d0  ;; go from 5 to 30 km/s

	;; array to save final values to
	vxcvals = dblarr(n_elements(linevelocities),n_elements(redshiftrange))

	;; start loop over linevelocities
	for kk=0,n_elements(linevelocities)-1 do begin

		vwidth = linevelocities[kk]
		result = templatespec(freqmin, freqmax, c_freq, vwidth)
		template_freq = reform(result[0,*])
		template_spec = reform(result[1,*])
		;; adjust the minimum
		template_spec = template_spec * 0.02
	
		crosscorrvals = dblarr(n_elements(redshiftrange))

		for jj=0,n_elements(redshiftrange)-1 do begin

			;; redshift the rest frame template to observed
			redshift = redshiftrange[jj]
			template_freq_obs = template_freq / (1d0 + redshift)


			sampled_template_spec = dblarr(n_elements(freqs))
			;; sample the template to have the same frequency coverage
			;;     remember that the channels are defined by the central frequency
			for ii=0,n_elements(freqs)-1 do begin

				chmin = freqs[ii] - chwidth/2d0
				chmax = freqs[ii] + chwidth/2d0
				freq_index = where(template_freq_obs ge chmin and template_freq_obs le chmax)
				;; find the contribution within a channel
				sampled_template_spec[ii] = total(template_spec[freq_index]) / n_elements(freq_index)
					
			endfor ;; index ii	

			;; Cross-correlate!!
			result = c_correlate(spec, sampled_template_spec, dblarr(n_elements(spec)))
			crosscorrvals[jj] = result[0]

			;; print sampled template to file
			ts = '/net/leubeek/data1/m82/crosscorr/templates/Ntemplate_spec_'+sigmastr+'_v'+strtrim(string(uint(vwidth)),2)+'.dat'
			openw,lun2,ts,/get_lun
			for ll=0,n_elements(freqs)-1 do printf,lun2,subbands[ll],' ',freqs[ll],' ',sampled_template_spec[ll],format='(D,A,D,A,D)'
			close,lun2
			free_lun,lun2
			

		endfor ;; index jj
		vxcvals[kk,*] = crosscorrvals

		plotfile = 'crosscorr/Ncrosscorr_' + sigmastr + '_v' + cgnumber_formatter(vwidth,decimal=0) + 'eps'
	        device,filename=plotfile,/encap
	        plot,redshiftrange,crosscorrvals,thick=3,charthick=3,xtitle='Redshift',ytitle='Cross-correlation',charsize=1,position=[0.15,0.15,0.9,0.9]
	        device,/close
		
	endfor ;; index kk

	;; print to file so you don't have to run it again
	logfile = 'Ncrosscorr_'+sigmastr+'.dat'
	openw,lun,logfile,/get_lun
	printf,lun,'redshift ',strtrim(string(uint(linevelocities)),1)
	for kk=0,n_elements(redshiftrange)-1 do begin
		ss = cgnumber_formatter(redshiftrange[kk],decimal=2)
		for jj=0,n_elements(linevelocities)-1 do begin
			ss = ss + ' ' + cgnumber_formatter(vxcvals[jj,kk],decimal=5)
		endfor ;; index jj
		printf,lun,ss
	endfor ;; index kk
	close,lun
	free_lun,lun
	

	plotfile = 'Ncrosscorr_'+sigmastr+'.eps'
	;; determine yrange
	ymin = min(vxcvals) - (0.1*min(vxcvals) + 0.1*max(vxcvals))/2d0
	ymax = max(vxcvals) + (0.1*min(vxcvals) + 0.1*max(vxcvals))/2d0
	yr = [ymin, ymax]
	sunset_colors,n_colors=254
	colors = intarr(n_elements(linevelocities))+1
	for kk=1,n_elements(colors)-1 do colors[kk] = colors[kk-1] + 253/(n_elements(colors)-1)

	device,filename=plotfile,/encap,/color
	plot, [-1,-1], [-1,-1],thick=4,ytitle='Cross-correlation',charsize=1,position=[0.15,0.15,0.8,0.9],charthick=3, color=0, xrange=[min(redshiftrange),max(redshiftrange)],yrange=yr
	for kk=0,n_elements(linevelocities)-1 do oplot, redshiftrange, reform(vxcvals[kk,*]), thick=4, color=colors[kk]
	names = strtrim(string(uint(linevelocities)),2)
	labelindex = [indgen(n_elements(names)/4)*4, n_elements(names)-1]
	cgcolorbar, divisions=n_elements(labelindex)-1, ticknames=names[labelindex], /vertical, position=[0.82,0.15,0.85,0.9],/right, ncolors=254,ticklen=0,bottom=1, charsize=0.8, charthick=3
	cblabel = textoidl('Line Width [km s^{-1}]')
	cgtext, cblabel, 0.9, 0.55, /normal, orientation=270, alignment=0.5, charsize=0.9, charthick=3
	cgtext, 'Redshift', 0.475, 0.05, /normal, alignment=0.5, charsize=1, charthick=3
	device,/close

	print, '======================================='
	print, '         redshift xc values            '
	print, '---------------------------------------'
	print, '  Line width    Redshift               '
	for kk=0,n_elements(linevelocities)-1 do print, '  '+cgnumber_formatter(linevelocities[kk],decimal=1)+' km/s     '+cgnumber_formatter(redshiftrange[where(reform(vxcvals[kk,*]) eq max(vxcvals[kk,*]))],decimal=2)

print,'done.'
end
