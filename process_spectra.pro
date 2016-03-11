pro process_spectra

;; start from new directory with only restored.corr images in it

	;; create fits files -- COMMENT ONE OF THESE OUT!
	spawn, 'casapy --nologger -c /home/morabito/scripts/RRLspectra/makefits.py' ;; no smoothing
	;spawn, 'casapy --nologger -c /home/morabito/scripts/RRLspectra/makefits_sm.py' ;; with smoothing

	;; clean up after casa
	if file_test('casapy*log') then file_delete, 'casapy*log'
	if file_test('ipython*log') then file_delete, 'ipython*log'
	if file_test('exportfits.last') then file_delete, 'exportfits.last'

	;; check if files are sorted and sort them if not
	if file_test('channel_images/SB*',/directory) ne 1 then begin
		;; first get a list of the files
		allimages = file_search('channel_images/*fits')
		sbs = []
		for ii=0,n_elements(allimages)-1 do begin
			tmp = strsplit(allimages[ii],'SB',/regex)
			sbs = [sbs, 'SB'+strmid(allimages[ii],tmp[1],3)]
		endfor
		sblist = sbs[rem_dup(sbs)]
		for ii=0,n_elements(sblist)-1 do begin
			sbdir = 'channel_images/'+sblist[ii]
			file_mkdir, sbdir
			file_move, 'channel_images/*'+sblist[ii]+'*fits', sbdir, /require_directory
		endfor
	endif

	;; run the extraction 
	runex = 0 ;; set to 0 if you want to run the extraction
	if runex eq 0 then begin
		print, 'STARTING SPECTRAL EXTRACTION ...'
		aperextraction
		print, 'SPECTRAL EXTRACTION COMPLETE.'
		print, '============================='
		print, ' '
	endif

	;; what is the smallest and largest beam size?
	readcol,'beamsizes.dat',sb,bmaj,bmin,bpa,format='A,D,D,D'
	beamsizes = sqrt(bmaj^2d0 + bmin^2d0)
	print, ' '
	print, '============================='
	print, ' '
	print, 'Plotting the Gaussian fits to identify bad subbands.'
	readcol,'aperextraction_fitparms.dat',subband,constant,height,widthx,widthy,centerx,centery,theta,format='A,D,D,D,D,D,D,D'
	pltsb = intarr(n_elements(subband))
	for ii=0,n_elements(pltsb)-1 do pltsb[ii] = fix(strmid(subband[ii],2))


	cgps_open,filename='aperextraction_fitparms.ps'
        cgloadct,0
        !p.multi=[0,3,3]
        plot,pltsb,constant,title='constant term',xtitle='subband',color=0
        plot,pltsb,height,title='scale factor',xtitle='subband'
        plot,pltsb,widthx,title='width x',xtitle='subband'
        plot,pltsb,widthy,title='width y',xtitle='subband'
        plot,pltsb,centerx,title='center x',xtitle='subband',yrange=[20,60]
        plot,pltsb,centery,title='center y',xtitle='subband',yrange=[20,60]
        ;; convert the radians to degrees
        thetadeg = theta * 360d0 / (2d0*!dpi)
        plot,pltsb,thetadeg,title='theta (deg)',xtitle='subband'
        cgps_close


	if file_test('sbplot',/directory) ne 1 then file_mkdir, 'sbplot'

	print, 'STARTING SUBBAND PLOTS ...'
	apersbplot
	print, 'SUBBAND PLOTS COMPLETE.'
	print, '============================='
	print, ' '
	print, 'PLOTTING CURVE-OF-GROWTH ...'
	curveofgrowth
	print, 'CURVE-OF-GROWTH COMPLETE.'
	print, '============================='
	print, ' '


	print, '============================='
	print, '    SUMMARY OF PROCESSING'
	print, '-----------------------------'
	print, 'Processed ',strtrim(fix(n_elements(rem_dup(sb))),2),' subbands'
	print, 'Aperture extraction files saved to aperextraction/'
	print, 'SB inspection plots saved to sbplot/'
	readcol,'imstddev.dat',sb,imstdev,format='A,D'
	failed = where(imstdev gt median(imstdev)*5)
	if total(failed) ge 0 then begin
		print, 'There were ',strtrim(fix(n_elements(failed)),2),' subbands that failed.'
	endif
	print, ' '
        print, 'Smallest beam size: '
        print, '-------------------'
        print, 'minor axis:',bmin[where(beamsizes eq min(beamsizes))]
        print, 'major axis:',bmaj[where(beamsizes eq min(beamsizes))]
        print, 'Largest beam size: '
        print, '-------------------'
        print, 'minor axis:',bmin[where(beamsizes eq max(beamsizes))]
        print, 'major axis:',bmaj[where(beamsizes eq max(beamsizes))]
	print, '============================='
	print, ' '
	print, 'Ready to run stack.pro'

print,'done.'
end
