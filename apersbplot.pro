pro apersbplot
;+
; NAME:
;     SBPLOT
; PURPOSE:
;     Plot the spectra extracted with gaussextraction.pro to plot,
;       calculate rms, etc. to filter bad subbands.
; EXPLANATION:
;     Reads the extracted spectra, fits a polynomial to the spectra,
;       calculates rms, median, etc. to identify bad subbands
;
; NOTES:
;       (1) Writes the following files:
;           - postscript file with spectra for each subband and poly fit
;           - .dat file that contains statistics on each subband
;
;       (2) Only writes spectra for good subbbands to .dat file. 
;               Filtering is done based on statistics calculated
;
;   MODIFICATION HISTORY:
;	Written in IDL by Leah K. Morabito October 2013
;	Revised 8 November 2013 (mostly comments)
;	Last revise 12 March 2014 (adapt for general pipeline use)
;-

	;; get a list of apertures
	aperfiles = file_search('aperextraction/*dat')

	for ii=0,n_elements(aperfiles)-1 do begin


		;;-- INPUTS
		readcol,aperfiles[ii],sb,freq,spec,format='A,D,D'
		sblist = sb[rem_dup(sb)]

		;;--SET UP OUTPUT FILES
		outdir = 'sbplot/'
		tmp = strsplit(aperfiles[ii],'/',/extract)
		plotfile = outdir + strmid(tmp[1],0,3) + '.ps'
		logfile = outdir + tmp[1]
		!p.thick=3.0
		!x.thick=3.0
		!y.thick=3.0
		!p.multi=[0,1,2]

		;;-- CREATE FINAL VECTORS

		nsubbands = n_elements(sblist)

        	clipchans = 5
		clipedge = clipchans - 3

		;;-- BEGIN LOOP OVER SUBBANDS

		;; start the log file
		openw,lun,logfile,/get_lun
		printf,lun,'sbid contlev rmsresid dynrange',format='(A)'

		;; start the plot file
		cgps_open,filename=plotfile,/landscape,/color
		cgloadct,39

		nsb = 0
		for jj=0,nsubbands-1 do begin

			;; locate the subband
			subband = sblist[jj]
			sbindex = where(sb eq subband)

			;; extract the spectral information
			freq_sb = freq[sbindex]/1d6  ;; convert to MHz
			spec_sb = spec[sbindex]

			;; clip the edges
			freq_sb = freq_sb[clipedge:n_elements(freq_sb)-clipedge]
			spec_sb = spec_sb[clipedge:n_elements(spec_sb)-clipedge]
		
			;;;;; SPEC
			;; linear fit
			coeff = linfit(freq_sb, spec_sb,yfit=yfit,sigma=stuff)

			;; calculate statistics and format for writing to file
			cont_level = median(spec_sb[clipedge:n_elements(spec_sb)-1-clipedge])
			rms_resid = stddev(spec_sb[clipedge:n_elements(spec_sb)-1-clipedge])

			dyn_range = rms_resid / cont_level
			cont_level = cgnumber_formatter(cont_level,decimal=5)
			rms_resid = cgnumber_formatter(rms_resid,decimal=5)
			dyn_range = cgnumber_formatter(dyn_range,decimal=5)
			printf,lun,subband,' ',cont_level,' ',rms_resid,' ',dyn_range,format='(A,A,A,A,A,A,A)'

			;; plot
			spec_fit = (spec_sb/yfit) -1d0
			yrange = [min(spec_sb)-1,max(spec_sb)+1]
			st = 'SB'+strtrim(string(subband),2)
			plot,freq_sb,spec_sb,psym=10,title=st,yrange=yrange,/ys
			oplot,freq_sb,yfit,linestyle=2,color=250
			plot,freq_sb,spec_fit,title=st,psym=10

		endfor ;index jj, loop over subbands

		close,lun
		free_lun,lun

		cgps_close

	endfor ;; index ii, loop over apertures

print, 'done.'
end
