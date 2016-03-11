pro aperextraction
;+
; NAME:
;     APEREXTRACTION
; PURPOSE:
;     For each subband, extract the spectra within a Gaussian area 
;	fitted to the source in the center.
;      
; EXPLANATION:
;     Reads fits cubes of data, fits a Gaussian to each subband,
;       extracts total spectra from the area.
;
; NOTES:
;       (1) Writes the following files:
;           - an .ps file with images of each channel and Gaussian regions
;           - a .dat file that contains a list of the extracted spectra
;	    - an .ps file with gaussian fit parameters for all subbands
;	    - a .dat file that contains a log of the gaussian fit parameters
;
;       (2) Only writes spectra for good subbbands to .dat file. 
;		Filtering is done based on successful Gaussian fits.
;
;   MODIFICATION HISTORY:
;      Written in IDL by Leah K. Morabito October 2013
;      Revised 8 November 2013 (mostly comments)
;      Revised 12 March 2014 (adapt for general pipeline use)
;      Last revised 14 January 2015 (minor cleanup of typos, etc.)
;-

	;;-- GET A LIST OF SUBBANDS
	sblist = file_search('channel_images/SB*')

        short = 0 ; >0 to test on only a few subbands, 0 to run the full extraction

        if short gt 0 then begin
		print, 'running test on only a few subbands and channels.'
                index = [0,4,59,65]
                sblist = sblist[index]
        endif

	;;-- OPEN THE GAUSSIAN FIT PARAMS LOGFILE
	openw,lun1,'aperextraction_fitparms.dat',/get_lun
        printf,lun1,'subband constant height widthx widthy centerx centery theta'

	if file_test('aperextraction',/directory) ne 1 then file_mkdir, 'aperextraction'

	;;sigmas = (dindgen(51)+10d0)/10d0
	;;sigmas = (dindgen(71)+10d0)/10d0
	sigmas = (dindgen(21)+10d0)/10d0
	for kk=0,n_elements(sigmas)-1 do begin
		sigfile = 'aperextraction/'+cgnumber_formatter(sigmas[kk],decimal=1)+'.dat'
		openw,lun2,sigfile,/get_lun
		printf,lun2,'subband frequency spectra'
		close,lun2
		free_lun,lun2
	endfor ;; index kk
	
	;;-- OPEN A LOGFILE TO SAVE BEAM SIZES
	openw,lun3,'beamsizes.dat',/get_lun
	printf,lun3,'subband bmaj bmin bpa'

	;;-- START LOOP OVER SUBBANDS
	sbs = []
	sbstddev = dblarr(n_elements(sblist))
	for jj=0,n_elements(sblist)-1 do begin

                tmp = strsplit(sblist[jj],'/',/extract)
                subband = tmp[1]
		sbs = [sbs, subband]

                print, 'Processing '+subband

                ;;-- get list of channel images and sort them properly
                chanlist = file_search('channel_images/'+subband+'/*fits')
                sortindex = intarr(n_elements(chanlist))
                for ii=0,n_elements(chanlist)-1 do begin
                        tmp = strsplit(chanlist[ii],'chan',/extract,/regex)
                        tmp1 = strsplit(tmp[1],'.',/extract)
                        sortindex[ii] = fix(tmp1[0])
                endfor ;; index ii
                chanlist = chanlist[sort(sortindex)]

                if short gt 0 then begin
                        chindex = [0,8,16,24]
                        chanlist = chanlist[chindex]
                endif


                ;;-- MAKE A CUBE
                fitsheader = headfits(chanlist[0])
		npixx = sxpar(fitsheader,'NAXIS1')
		npixy = sxpar(fitsheader,'NAXIS2')
                xcen = sxpar(fitsheader,'NAXIS1') / 2d0
                ycen = sxpar(fitsheader,'NAXIS2') / 2d0
                xmin = xcen - 50 
                xmax = xcen + 50
                ymin = ycen - 50
                ymax = ycen + 50
                datacube = dblarr(n_elements(chanlist),npixx,npixy)
                for ii=0,n_elements(chanlist)-1 do datacube[ii,*,*] = mrdfits(chanlist[ii],0,h,/silent)

		;;-- squash the cube and cut off the noisy edges
		;;-- faster to do this than to create a mask
		tmp = total(datacube,1)
		contimage = tmp[xmin:xmax,ymin:ymax]

		sbstddev[jj] = stddev(tmp)
				
		;;-- fit a gaussian to define an aperture
                g2dparms1 = [0d0, max(contimage), 3, 3, xcen-xmin, ycen-ymin, 0d0]
                parinfo = replicate({value:0.D,fixed:0, limited:[0,0], limits:[0.D,0]},7)
                parinfo[*].value = g2dparms1
                parinfo[*].fixed = 0
                result = mpfit2dpeak(contimage, a, parinfo=parinfo, estimates=g2dparms1, perror=perr, bestnorm=bn,/tilt)
;		result = gauss2dfit(contimage, a, /tilt)
		;;-- normalize
                result = result - min(result)
		normalizeresult = result / max(result)

		;;-- use sigma cutoffs to make masks for apertures of varying sizes
		sbmasks = dblarr(n_elements(sigmas),npixx,npixy)
		for ii=0,n_elements(sigmas)-1 do begin

			sigma = sigmas[ii]
			tmp = dblarr(size(normalizeresult,/dimensions))
			tmp[where(normalizeresult ge (1d0-erf(sigma/sqrt(2d0))))] = 1d0
			;; zero pad
			sbmasks[ii,xmin:xmax,ymin:ymax] = tmp
		endfor ;; index jj, loop over sigmas

		plotfile = 'aperextraction/'+subband+'.ps'
		cgps_open,filename=plotfile,/color
	        !P.charsize=1.5
		loadct,5
		
		;;-- NOW LOOP OVER CHANNELS TO:
			;; (b) use an aperture to extract the spectra
		stuff = 1
		if stuff gt 0 then begin
		for ii=0,n_elements(chanlist)-1 do begin

			;;-- GET BEAM INFORMATION
			;;-- open the fits header 
                        fitsheader = headfits(chanlist[ii])
                        freq = sxpar(fitsheader,'CRVAL4')
                        pixscale = abs(sxpar(fitsheader,'CDELT1'))
                        ;;-- beam size, convert to arcsec
                        bmaj = sxpar(fitsheader,'BMAJ') * 60d0 * 60d0 
                        bmin = sxpar(fitsheader,'BMIN') * 60d0 * 60d0
			;;-- beam position angle
                        bpa = sxpar(fitsheader,'BPA')  ;; degrees
			printf,lun3,subband,' ',cgnumber_formatter(bmaj,decimal=4),' ',cgnumber_formatter(bmin,decimal=4),' ',cgnumber_formatter(bpa,decimal=4),format='(7A)'

                        chanimage = reform(datacube[ii,xmin:xmax,ymin:ymax])

			;; extract from apertures 
			!p.multi=[0,8,8]
			for kk=0,n_elements(sigmas)-1 do begin

	                        ;;-- plotting
        	                dummy = [0,0]
                	        thispos = [0.1,0.1,0.9,0.9]

				;;-- blank plo to label the subband
	                        if kk eq 0 then begin
        	                        cgloadct, 0
                	                plot,dummy, xstyle=4, ystyle=4, color=255,Position=thispos
                        	        sbtitle = 'ch '+strtrim(string(ii),2)
                                	cgText,0.03,0.93,sbtitle,/normal,color=0,charsize=0.8,charthick=3
	                                cgloadct, 33
					cgimage, chanimage
        	                endif

		                sigfile = 'aperextraction/'+cgnumber_formatter(sigmas[kk],decimal=1)+'.dat'

				mask = reform(sbmasks[kk,xmin:xmax,ymin:ymax])

				tmp = chanimage * mask
				spec = total(tmp)	

	                        ;;-- PLOT THE CHANNEL
                                cgimage, chanimage 
                                cgloadct, 0
                                cgcontour, mask, /onimage, levels=[0,1], color=255,label=0
                                cgloadct, 33

				openw,lun2,sigfile,/append,/get_lun
				printf,lun2,subband,' ',cgnumber_formatter(freq,decimal=4),' ',cgnumber_formatter(spec,decimal=7),format='(A,A,A,A,A)'
				close,lun2
				free_lun,lun2

			endfor ;; index kk, loop over sigmas

		endfor ;; index ii, loop over channels
		endif
		cgps_close

		printf,lun1,subband,' ',a[0],' ',a[1],' ',a[2],' ',a[3],' ',a[4],' ',a[5],' ',a[6],format='(A,A,D,A,D,A,D,A,D,A,D,A,D,A,D)'

	endfor ;; index jj, loop over subbands

        close,lun1
	free_lun,lun1

	close,lun3
	free_lun,lun3
	
	openw,lun1,'imstddev.dat',/get_lun
	for jj=0,n_elements(sbstddev)-1 do printf,lun1,sbs[jj],' ',sbstddev[jj],format='(A,A,D)'
	close,lun1
	free_lun,lun1

print, 'done.'
end
