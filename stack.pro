pro stack, sigma, redshift
;+
; NAME:
;     STACK
; PURPOSE:
;     Find subbands with CRRLs and stack them.
; EXPLANATION:
;     Reads the extracted spectra, filters based on results of 
;	sbplot.pro, fits the bandpass of each SB (with an option
;	to use previous fits), converts to velocity space, and 
;	plots final stack.
;
; NOTES:
;       (1) Writes the following files:
;           - .ps file with pre and post bandpass subtraction, and 
;		stacked spectrum at the end
;           - .log file that contains information on each subband, 
;		including fit order for bandpass
;	    - .log file that contains the sorted, stacked spectrum only
;	    - .log file that contains the total spectrum with more information
;
;   MODIFICATION HISTORY:
;      Written in IDL by Leah K. Morabito October 2013
;      Revised 11 November 2013 to inspect the rms
;      Last revised 12 March 2014 (adapt for general pipeline use)
;-

	sigmastr = cgnumber_formatter(sigma,decimal=1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                   CONSTANTS                  ;;

        speedoflight = 2.99792458d5 ;; km/s
	sysvel = redshift * speedoflight ;; km/s

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;          BARYCENTRIC CORRECTION              ;;

	;; get date and time (UTC) of observation
	cubelist = file_search('channel_images/SB*/*fits')
	testcube = mrdfits(cubelist[0],0,h)
	utctime = sxpar( h, 'DATE-OBS')  
	radeg = sxpar(h, 'OBSRA')
	decdeg = sxpar(h, 'OBSDEC')
	tmp = strsplit(utctime,/extract,'-')
	yr = tmp[0]
	mo = tmp[1]
	tmp2 = strsplit(tmp[2],/extract,'T')
	day = tmp2[0]
	tmp3 = strsplit(tmp2[1],/extract,':')
	hr = tmp3[0]
	min = tmp3[1]
	sec = tmp3[2]
	jd = julday(mo,day,yr,hr,min,sec)
	baryvel, jd, 0, dvelh, dvelb, /jpl
	ra = radeg/!RADEG
	dec = decdeg/!RADEG
	vshift = dvelb[0]*cos(dec)*cos(ra) + dvelb[1]*cos(dec)*sin(ra) + dvelb[2]*sin(dec)
	print, 'Barycentric correction: ',cgnumber_formatter(vshift, decimal=3)



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;		build list of subbands		;;

	readcol,'aperextraction/' + sigmastr + '.dat',sb,freq,spec,format='A,D,D'
	sblist = sb[rem_dup(sb)]
	readcol,'sbplot/' + sigmastr + '.dat',sbid,contlev,rmsresid,dynrange,format='A,D,D,D'
	;; identify subbands based on rmsresid and dynrange that should not be used
	rmsmed = median(rmsresid)
	dynmed = median(dynrange)
	badsb = where(rmsresid gt 10*rmsmed or dynrange gt 10*dynmed, complement=goodsb)

	;; get rid of them from sblist
	for i=0,n_elements(badsb)-1 do begin
		tmpind = where(sblist eq sbid[badsb[i]],complement=goodtmpind)
		sblist = sblist[goodtmpind]
                print, 'Subband ',string(sbid[badsb[i]]),' was removed due to rms noise.'	
	endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                  list of rrls                ;;

	;; frequencies in MHz
        ;; ALPHA TRANSITIONS
        readcol, '/net/leubeek/data2/LLINE/RRL_CI_alpha.txt', ca1, ca2, cann, calpha_freq, format='A,A,D,D', skipline=1, /silent
        ;; BETA TRANSITIONS
        readcol, '/net/leubeek/data2/LLINE/RRL_CI_beta.txt', cb1, cb2, cbnn, cbeta_freq, format='A,A,D,D', skipline=1, /silent
        ;; GAMMA TRANSITIONS
        readcol, '/net/leubeek/data2/LLINE/RRL_CI_gamma.txt', cg1, cg2, cgnn, cgamma_freq, format='A,A,D,D', skipline=1, /silent


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;		SET UP OUTPUT FILES		;;

	if file_test('stack') ne 1 then file_mkdir, 'stack'
	outname = 'stack/' + sigmastr + '_z' + cgnumber_formatter(redshift,decimal=2)
	plotfile = outname + '.ps'
	logfile = outname + '.dat'
	specfile = outname + '.spec.dat'
	totalspecfile = outname + '.total.spec.dat'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;         FIT BANDPASS INTERACTIVELY?          ;;
;;           (must be done 1st time)		;;

        bpyesno=''
        read,bpyesno,prompt='y = fit the bandpass, n = use previous fits: '
        if bpyesno eq 'n' then begin
                readcol,logfile,fitsb,contlev,rmsresid,dynrange,polyord,format='A,D,D,D,D'
        endif else begin
		openw,lun,logfile,/get_lun
		printf,lun,'sbid contlev rmsresid dynrange polyord',format='(A)'
	endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;          SET UP PLOTTING PARAMETERS          ;;

	!p.thick=3.0
	!x.thick=3.0
	!y.thick=3.0
	!p.multi=[0,1,2]
	!P.font=0
	!x.margin=[12,8]
        dumx = [0,0]
        dumy = [0,0]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;              CREATE FINAL VECTORS            ;;

	nsubbands = n_elements(sblist)

        clipchans = 6
        clipedge = clipchans - 3


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;	    BEGIN LOOP OVER SUBBANDS		;;


	ps_start,filename=plotfile,/landscape,/color
	cgloadct,39

	stack_vel = 0d0
	stack_spec = 0d0
	stack_freq = 0d0
	stack_weights = 0d0
	alphafreqs = 0d0
	nsb = 0	
	sbnames = 0d0	

	for i=0,n_elements(sblist)-1 do begin

		subband = sblist[i]

		sbindex = where(sb eq subband)

		freq_sb = freq[sbindex]/1d6  ;; convert to MHz
		spec_sb = spec[sbindex]

		;; clip the edges
	;; 	this would be the corrrect way to clip
	;;	freq_sb_clip = freq_sb[clipedge:n_elements(freq_sb)-1-clipedge]
	;;	spec_sb_clip = spec_sb[clipedge:n_elements(spec_sb)-1-clipedge]

		
	;; 	but there is a missing channel (channel 32)
		freq_sb_clip = freq_sb[clipedge:n_elements(freq_sb)-clipedge]
		spec_sb_clip = spec_sb[clipedge:n_elements(spec_sb)-clipedge]

		;; test whether there's a line in the subband
		rest_freqs = freq_sb[clipchans:n_elements(freq_sb) - 1 - clipchans] * (1d0 + redshift)
	;;	rest_freqs = freq_sb[clipchans:n_elements(freq_sb) - clipchans] * (1d0 + redshift)
		alphafreqind = where(calpha_freq gt min(rest_freqs) and calpha_freq lt max(rest_freqs))
		tmp = calpha_freq[alphafreqind]
		alphafreq = tmp[0]
	
		;; start the loop if there's a line
		if total(alphafreqind) gt 0 then begin
			
			;; shift frequency to velocity
			velocity = ((alphafreq - freq_sb_clip) / freq_sb_clip * speedoflight) - sysvel + vshift

			;; blank line before fitting
			blankindex = where(velocity gt -100d0 and velocity lt 100d0) ;; assume line ~100 kms/s
			spec_blank = spec_sb_clip
			spec_blank[blankindex] = !values.f_nan
			spec_mean = mean(spec_blank[where(finite(spec_blank) ne 0)])
			spec_median = median(spec_blank[where(finite(spec_blank) ne 0)])
			spec_blank[blankindex] = spec_median
			
			;; FIT THE BANDPASS
			if bpyesno eq 'y' then begin

				;; plot for user to see
	                        set_plot,'x'
				!p.multi=0
	       	                plot, velocity, spec_blank, yrange=[min(spec_blank)-0.1d0*min(spec_blank), max(spec_blank)+0.1d0*max(spec_blank)], /ys
	              		tryfit: read,prompt='Please enter the order of the polynomial you want to try: ',order1
	                        ord = double(order1)
	               	        coeff = poly_fit(velocity, spec_blank, ord, /double)
	               	        ;coeff = robust_poly_fit(velocity, spec_blank, ord, /double)
                        	bp_poly = coeff[0]
	                        for k=1,n_elements(coeff)-1 do bp_poly = bp_poly + (coeff[k] * velocity^double(k) )
	
				;; overplot the fit for user to see
	       	                oplot, velocity, bp_poly, linestyle=2
                	        yesno = ''
	                       	read,prompt='Do you want to refit? (y/n): ', yesno
	                        if yesno eq 'y' then goto, tryfit

				;; once fitting is finished
				spec_bp = (spec_sb_clip / bp_poly) - 1d0

	                        ;; print statistics to logfile
	                        cont_level = median(spec_bp)
        	                rms_resid = stddev(spec_bp)
                	        dyn_range = rms_resid / cont_level
        		        weights = dblarr(n_elements(spec_bp))
                	        weights[*] = 1d0/rms_resid
        	                cont_level = cgnumber_formatter(cont_level,decimal=5)
                	        rms_resid = cgnumber_formatter(rms_resid,decimal=5)
                        	dyn_range = cgnumber_formatter(dyn_range,decimal=5)
	                        ord = cgnumber_formatter(ord,decimal=0)
        	                printf,lun,subband,' ',cont_level,' ',rms_resid,' ',dyn_range,' ',ord,format='(A,A,A,A,A,A,A,A,A)'

			endif else begin
				
				;; if using previous fits
				ord = polyord[where(sbid eq subband)]
				ord1 = ord[0]
				coeff = robust_poly_fit(velocity, spec_blank, ord1)
                                bp_poly = coeff[0]
                                for k=1,n_elements(coeff)-1 do bp_poly = bp_poly + (coeff[k] * velocity^double(k) )
	                        spec_bp = (spec_sb_clip / bp_poly) - 1d0
	
				rms_resid = stddev(spec_bp)
				wgt = 1d0/rms_resid
				
                                weights = replicate(wgt,n_elements(spec_bp))

			endelse
			
			;; plot to file	
                        set_plot,'ps'
                        !p.multi=[0,1,2]
			vt = 'Velocity [km/s]'
			nt = 'Normalized flux'

			yrange = [min(spec_sb_clip)-0.01*min(spec_sb_clip),max(spec_sb_clip)+0.01*max(spec_sb_clip)]
			pretitle = 'Subband '+string(subband)
			plot,velocity,spec_blank,title=pretitle,yrange=yrange,xtitle=vt,ytitle=nt,charsize=1
	;		plot,velocity,spec_sb_clip,title=pretitle,yrange=yrange,xtitle=vt,ytitle=nt,charsize=1
			oplot,velocity,bp_poly,linestyle=2,color=250
                        yrange = [min(spec_bp)-0.01,max(spec_bp)+0.01]
                        posttitle = 'Corrected for bandpass'
			plotcontinuum = bp_poly / bp_poly - 1d0
			result = linfit(velocity, spec_bp, yfit=yfit)
			plot,velocity,spec_bp,title=posttitle,yrange=yrange,linestyle=0,xtitle=vt,ytitle=nt,charsize=0.95
			oplot,velocity,yfit,linestyle=2,color=250

			stack_vel = [stack_vel, !values.f_nan, velocity]
			stack_spec = [stack_spec, !values.f_nan, spec_bp]
			stack_freq = [stack_freq, !values.f_nan, freq_sb_clip]
			stack_weights = [stack_weights, !values.f_nan, weights]

                        ;; book keeping
                        alphafreqs = [alphafreqs, !values.f_nan, replicate(alphafreq,n_elements(velocity))]
                        nsb = nsb + 1
			sbnum = strmid(subband,2,4)
			print, sbnum
                        sbnames = [sbnames, !values.f_nan, replicate(fix(sbnum),n_elements(velocity))]

						
		endif

	endfor ;; index i, loop over subbands

	print, strtrim(string(nsb),2), ' subbands had lines.'
	if bpyesno eq 'y' then begin
		close,lun
		free_lun,lun
	endif

	stack_vel = stack_vel[1:*]
	stack_spec = stack_spec[1:*]
	stack_freq = stack_freq[1:*]
	stack_weights = stack_weights[1:*]

	allvel = stack_vel
	allspec = stack_spec
	allfreq = stack_freq
	allweights = stack_weights
	allsb = sbnames[1:*]
	allalpha = alphafreqs[1:*]

	;; write the final TOTAL spectrum to file
	openw,lun1,totalspecfile,/get_lun
	printf,lun1,'Frequency  Spectra Weights Subband Alphafreq',format='(A)'
			
	for kk=0,n_elements(allspec)-1 do printf, lun1, allfreq[kk],' ', allspec[kk],' ',allweights[kk],' ',allsb[kk],' ',allalpha[kk],format='(D,A,D,A,D,A,D,A,D)'
	close,lun1
	free_lun,lun1

	;; clear out the NaNs
	stack_vel = stack_vel[where(finite(stack_vel) ne 0)]
	stack_spec = stack_spec[where(finite(stack_spec) ne 0)]

	stack_freq = stack_freq[where(finite(stack_freq) ne 0)]
	stack_weights = stack_weights[where(finite(stack_weights) ne 0)]
	stack_sb = allsb[where(finite(allsb) ne 0)]

	;; sort by velocity
	sortedvel = stack_vel[sort(stack_vel)]
        sortedspec = stack_spec[sort(stack_vel)]
        sortedfreq = stack_freq[sort(stack_vel)]
        sortedweights = stack_weights[sort(stack_vel)]
	sortedsb = stack_sb[sort(stack_vel)]

	;; calculate coverage
	indmin = n_elements(sortedvel)/2 - 50
	indmax = n_elements(sortedvel)/2 + 50
	fitvel = sortedvel[indmin:indmax]
        velx = dindgen(n_elements(sortedvel))
	fitvelx = velx[indmin:indmax]
	coeffs = linfit(fitvelx, fitvel, yfit=velvals)
	linvel = coeffs[0] + coeffs[1]*velx
	resid = sortedvel - linvel
	fitvel = coeffs[0] + coeffs[1]*fitvelx
	goodresid = where(resid ge -5d0 and resid le 5d0)
	velrange = sortedvel[min(goodresid):max(goodresid)]
	print, 'Use velocity range: ',min(velrange),' to ',max(velrange)
	xr = [-400,400]
	yr = [-0.006, 0.004]

	;; Plotting
        contline = dblarr(n_elements(sortedvel))
	goodvcov = where(sortedvel ge min(velrange) and sortedvel le max(velrange),complement=zerocov)
	
	;; Plot of total spectrum
	plot, sortedvel,sortedspec,linestyle=0,title='Stacked spectrum: Original Resolution',xtitle=vt
	oplot, sortedvel, contline, linestyle=2, color=250
	sortedspeccov = sortedspec
	sortedspeccov[zerocov] = 0d0
	;; Plot of specturm with good velocity coverage
	plot, sortedvel,sortedspeccov,linestyle=0,title='Stacked spectrum: Original Resolution, Good velocity coverage', xtitle=vt, xrange=xr
        oplot, sortedvel, contline, linestyle=2, color=250
	;; Smooth with boxcar and plot
	smoothed = smooth(sortedspec,35,/edge_truncate)
	smoothed[zerocov] = 0d0
	plot, sortedvel, smoothed, title='Moving boxcar, window 35', xrange=xr,yrange=yr
        oplot, sortedvel, contline, linestyle=2, color=250

	;; Use Savitzky-Golay order 1 filter
	savgolfilter = savgol(15, 15, 0, 1)
	savgolplot = convol(sortedspec, savgolfilter)
	savgolplot[zerocov] = 0d0
	;; weighted
	savgolwplot = convol(sortedspec * (sortedweights)/mean(sortedweights), savgolfilter)
	savgolwplot[zerocov] = 0d0
	plot, sortedvel, savgolwplot, title='S-G Filter, order 1, window 31, weighted', xrange=xr,yrange=yr
        oplot, sortedvel, contline, linestyle=2, color=250

	;; Use Savitzky-Golay order 2 filter
	savgolfilter2 = savgol(26, 26, 0, 2)
	savgolwplot2 = convol(sortedspec*sortedweights / mean(sortedweights), savgolfilter2)
	savgolwplot2[zerocov] = 0d0
        plot, sortedvel, savgolwplot2, title='S-G Filter, order 2, window 53, weighted', xrange=xr,yrange=yr
        oplot, sortedvel, contline, linestyle=2, color=250
	gls = lowess(sortedvel, sortedspec, 53, 1, noise)
	gls[zerocov] = 0d0

	;; General Least-Squares filter
	plot, sortedvel, gls, title='General Least Squares filter, window 53', xrange=xr,yrange=yr
	oplot, sortedvel, contline, linestyle=3, color=250


	;; Plot velocity coverage information
	!P.multi=0
	toppos = [0.2,0.6,0.9,0.9]
	midpos = [0.2,0.4,0.9,0.6]
	botpos = [0.2,0.15,0.9,0.4]
	xr=[-300,300]
        plot, sortedvel, savgolwplot, title='S-G Filter, window 31, weighted by 1/SD', xrange=xr,yrange=yr,position=toppos,XTickformat='(A1)'
        oplot, sortedvel, contline, linestyle=2, color=250
	dv = textoidl(" \Delta v")
	velcov = 0 -smooth(sortedvel[1:*] - sortedvel[0:n_elements(sortedvel)-2], 31)
	velcov = [velcov[0], velcov]
	plot, sortedvel, velcov, ytitle=dv, xrange=xr,position=midpos,/noerase, XTickformat='(A1)'
	coeffs = poly_fit(sortedvel, velcov, 12)
	velcovfit = coeffs[0]
	for k=1,n_elements(coeffs)-1 do velcovfit = velcovfit + (coeffs[k] * sortedvel^double(k) )
	oplot,sortedvel, velcovfit, color=250,thick=4

	plothist, sortedvel, /autobin, peak=1, position=botpos, /noerase, xtitle=vt, ytitle='normalized density', xrange=xr, yrange=[0.4,1.05],ystyle=1,thick=4
	ninetyfivepc = [0.95,0.95]
	oplot, xr, ninetyfivepc, linestyle=1,color=250,thick=4

	velones = dblarr(n_elements(sortedvel))
	velones[*] = 1d0

	!P.multi=[0,1,2]	
	velx = dindgen(n_elements(sortedvel))
	plot,dumx,dumy,xrange=[min(velx),max(velx)],yrange=[min(sortedvel),max(sortedvel)],ytitle=vt,title='Velocity coverage'
	oplot,velx,sortedvel,thick=4,color=250
	oplot,velx,linvel,linestyle=2,color=50,thick=4

	plot,velx,resid,title='residuals',color=0,thick=4
	oplot,velx[goodresid],resid[goodresid],color=250,thick=5

	;; end plotting!
	ps_end

	;; write spectra to file
	openw,lun2,specfile,/get_lun
	printf,lun2,'velocity  unsmoothed  smoothed35 sgsord1 sgord1w sgord2 sgord2w weights',format='(A)'
        smoothed = smooth(sortedspec,35)
	sgord1 = convol(sortedspec, savgolfilter)
	sgord1w = convol(sortedspec*sortedweights/mean(sortedweights), savgolfilter)
	sgord2 = convol(sortedspec, savgolfilter2)
	sgord2w = convol(sortedspec*sortedweights/mean(sortedweights), savgolfilter2)
	for i=0,n_elements(sortedvel)-1 do printf,lun2,sortedvel[i],' ',sortedspec[i],' ',smoothed[i],' ',sgord1[i],' ',sgord1w[i],' ',sgord2[i],' ',sgord2w[i],' ',sortedweights[i],format='(D,A,D,A,D,A,D,A,D,A,D,A,D,A,D)'
	close,lun2
	free_lun,lun2


	print,'Processed a total of ',nsubbands,'.'
	print,'There were ',nsb,' subbands stacked.'
	finalnoise = 1d0 / sqrt(total(sortedweights[where(sortedvel gt -150 and sortedvel lt 150)]))
	print,'The noise in the final spectrum is: ', cgnumber_formatter(finalnoise,decimal=3)
	linepeak = min(sortedspec[where(sortedvel gt -50 and sortedvel lt 50)])
	print,'The raw spectrum peaks at ', cgnumber_formatter(linepeak,decimal=3)
	print,'For an SNR of: ',cgnumber_formatter(abs(linepeak/finalnoise),decimal=3)

	set_plot,'x'
stop

print, 'done.'
end
