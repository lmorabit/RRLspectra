pro newcurveofgrowth

        extractionfiles = file_search('aperextraction/*dat')
        tmpfile = extractionfiles[n_elements(extractionfiles)-1]
        tmp = strsplit(tmpfile,'/',/extract)
        maxaper = strmid(tmp[1],0,3)
        maxaper = double(maxaper)
	maxaperarea = !dpi*maxaper^2d0

        napers = n_elements(extractionfiles)

        aper = dblarr(napers)
        spectotal = dblarr(napers)
	aperarea = dblarr(napers)

        for j=0,napers-1 do begin

                exfile = extractionfiles[j]
                tmp = strsplit(exfile,'/',/extract)
                aper[j] = strmid(tmp[1],0,3)
                readcol,exfile,subband,freq,spec,format='A,D,D',/silent

                sblist = subband(rem_dup(subband))

                tmpfreq = 0d0
                tmpspec = 0d0
                for i=0,n_elements(sblist)-1 do begin
                        sbindex = where(subband eq sblist[i])
                        tmpfreq = [tmpfreq, freq[sbindex]]
                        tmpspec = [tmpspec, total(spec[sbindex])/n_elements(sbindex)]
                endfor

                freq = tmpfreq[1:n_elements(tmpfreq)-1]
                spec = tmpspec[1:n_elements(tmpspec)-1]

                nsb = n_elements(sblist)
                spectotal[j] = total(spec)/nsb ;;* ( (!dpi*aper[j]^2d0) / minaperarea)
		aperarea[j] = !dpi*aper[j]^2d0

        endfor ;index j
        ;spectotal = spectotal / max(spectotal) ;; normalize
        ;spectotal = spectotal / spectotal[14] ;; normalize

	stop

	!P.multi=0
        pltfl = 'curve_of_growth.ps'
        cgps_open, filename=pltfl, /color
        xt = textoidl("Aperture Radius")
        yt = 'Arbitrary Units'
        plot,aper,spectotal,charsize=2, xrange=[1.8,max(aper)-0.5],/xs, /ys, xtitle=xt, ytitle=yt
        cgps_close

end


