t#pragma rtGlobals=1		// Use modern global access method.

function trackpeak(wpre,wlen,outpeak,start,stop,imin,imax)
	string wpre
	wave outpeak,wlen
	variable start, stop,imax,imin
	variable i
	wave W_coef
	wave W_sigma
	string tempwavey
	for (i = imin; i < imax;i += 1)
		tempwavey=wpre+num2str(i)
		wavestats/Q/R=[start,stop] $tempwavey
		outpeak[i-1]=wlen[V_maxloc]		
	endfor
end

function trackgauss(wpre,outpeak,outerror,start,stop,width,imax)
	string wpre
	wave outpeak,outerror
	variable start, stop, width,imax
	
	variable i
	wave W_coef
	wave W_sigma
	string tempwavey, tempwavex
	tempwavex=wpre+"0"
	wave tempx=$tempwavex
	for (i = 0; i < imax;i += 1)
		tempwavey=wpre+num2str(i)
		wavestats/Q/R=[start,stop] $tempwavey
		CurveFit/Q/NTHR=0 gauss,  $tempwavey [V_maxloc-width/2, V_maxloc+width/2] /X=$tempwavex /D 
		outpeak[i-1]=W_coef[2]
		outerror[i-1]=W_sigma[2]		
	endfor
end

function trackgauss2(wpre,stnum,wlen,outpeak,outerror,start,stop,width)
	string wpre
	wave outpeak,outerror,wlen
	variable stnum,start, stop, width
	
	variable i,imax
	wave W_coef
	wave W_sigma
	string tempwavey
	wavestats/Q outpeak
	imax=V_npnts+stnum
	for (i = stnum; i < imax;i += 1)
		tempwavey=wpre+num2str(i)
		wavestats/Q/R=[start,stop] $tempwavey
		CurveFit/Q/NTHR=0 gauss,  $tempwavey [V_maxloc-width/2, V_maxloc+width/2] /X=wlen /D 
		outpeak[i-stnum]=W_coef[2]
		outerror[i-stnum]=W_sigma[2]		
	endfor
end

function trackgauss3(wpre, stnum, increment, endnum, wlen,suffix, temperature,timing,start,stop,width)
	string wpre,suffix
	wave temperature, timing,wlen
	variable stnum,endnum,increment, start, stop, width
	variable i,imax
	string tempwavey,peakname,errname,tempname,timename
	imax=floor((endnum-stnum)/increment)+1
	peakname="peak"+suffix
	errname="err"+suffix
	tempname="temp"+suffix
	timename="time"+suffix
	Make/N=(imax) /D $peakname
	Make/N=(imax) /D $errname
	Make/N=(imax) /D $tempname
	Make/N=(imax) /D $timename
	wave peakuse=$peakname
	wave erruse=$errname
	wave tempuse=$tempname
	wave timeuse=$timename
	for (i = 0; i < imax;i += 1)
		tempwavey=wpre+num2str(i*increment+stnum)
		wavestats/Q/R=[start,stop] $tempwavey
		CurveFit/Q/N=1/NTHR=0 gauss,  $tempwavey [V_maxloc-width/2, V_maxloc+width/2] /X=wlen /D 
		wave wcf=W_coef
		wave wsg=W_sigma 
		peakuse[i]=wcf[2]
		erruse[i]=wsg[2]
		tempuse[i]=temperature[i*increment+stnum-1]
		timeuse[i]=timing[i*increment+stnum-1]	
		killwaves wcf,wsg
	endfor
	//timeuse=timeuse/60
	display peakuse vs timeuse
	appendtograph/R /C=(0,0,65535) tempuse vs timeuse
	Label left "Peak Wavelength (nm)";DelayUpdate
	Label bottom "Time, seconds";DelayUpdate
	Label right "Temperature, ¼C";DelayUpdate
	Legend/C/N=text0/F=0/H={30,4,10}/A=MC
	ModifyGraph grid(left)=1,grid(bottom)=1
	ModifyGraph mirror(bottom)=2
end

function trackgauss3min(wpre, stnum, increment, endnum, wlen,suffix, temperature,timing,start,stop,width)
	string wpre,suffix
	wave temperature, timing,wlen
	variable stnum,endnum,increment, start, stop, width
	variable i,imax
	string tempwavey,peakname,errname,tempname,timename
	imax=floor((endnum-stnum)/increment)+1
	peakname="peak"+suffix
	errname="err"+suffix
	tempname="temp"+suffix
	timename="time"+suffix
	Make/N=(imax) /D $peakname
	Make/N=(imax) /D $errname
	Make/N=(imax) /D $tempname
	Make/N=(imax) /D $timename
	wave peakuse=$peakname
	wave erruse=$errname
	wave tempuse=$tempname
	wave timeuse=$timename
	for (i = 0; i < imax;i += 1)
		tempwavey=wpre+num2str(i*increment+stnum)
		wave tempwaveyuse=$tempwavey
		duplicate tempwaveyuse tmpneg
		tmpneg=-1*tempwaveyuse
		wavestats/Q/R=[start,stop] tmpneg
		CurveFit/Q/N=1/NTHR=0 gauss,  tmpneg [V_maxloc-width/2, V_maxloc+width/2] /X=wlen /D 
		wave wcf=W_coef
		wave wsg=W_sigma 
		peakuse[i]=wcf[2]
		erruse[i]=wsg[2]
		tempuse[i]=temperature[i*increment+stnum-1]
		timeuse[i]=timing[i*increment+stnum-1]	
		killwaves wcf,wsg
	endfor
	//timeuse=timeuse/60
	display peakuse vs timeuse
	appendtograph/R /C=(0,0,65535) tempuse vs timeuse
	Label left "Peak Wavelength (nm)";DelayUpdate
	Label bottom "Time, seconds";DelayUpdate
	Label right "Temperature, ¼C";DelayUpdate
	Legend/C/N=text0/F=0/H={30,4,10}/A=MC
	ModifyGraph grid(left)=1,grid(bottom)=1
	ModifyGraph mirror(bottom)=2
end

function trackgauss4(wpre, stnum, increment, endnum, wlen,suffix, temperature,timing,start,stop,width) //with normalization
	string wpre,suffix
	wave temperature, timing,wlen
	variable stnum,endnum,increment, start, stop, width
	variable i,imax,maxloc,leftminloc,rightminloc,yleft,yright,slope,intcpt
	string tempwavey,peakname,errname,tempname,timename
	imax=floor((endnum-stnum)/increment)+1
	peakname="peak"+suffix
	errname="err"+suffix
	tempname="temp"+suffix
	timename="time"+suffix
	Make/N=(imax) /D $peakname
	Make/N=(imax) /D $errname
	Make/N=(imax) /D $tempname
	Make/N=(imax) /D $timename
	wave peakuse=$peakname
	wave erruse=$errname
	wave tempuse=$tempname
	wave timeuse=$timename
	for (i = 0; i < imax;i += 1)
		tempwavey=wpre+num2str(i*increment+stnum)
		wavestats/Q/R=[start,stop] $tempwavey
		maxloc=V_maxloc
		wavestats/Q/R=[maxloc-(1.5*width),maxloc-(0.5*width)] $tempwavey
		leftminloc=V_minloc
		yleft=V_min
		wavestats/Q/R=[maxloc+(0.5*width),maxloc+(1.5*width)] $tempwavey
		rightminloc=V_minloc
		yright=V_min
		slope=(yright-yleft)/(rightminloc-leftminloc)
		intcpt=yright-slope*rightminloc
		duplicate $tempwavey tempwaveynorm
		tempwaveynorm=tempwaveynorm/(slope*x+intcpt)
		CurveFit/Q/N=1/NTHR=0 gauss,  tempwaveynorm [maxloc-width/2, maxloc+width/2] /X=wlen /D
		wave wcf=W_coef
		wave wsg=W_sigma 
		peakuse[i]=wcf[2]
		erruse[i]=wsg[2]
		tempuse[i]=temperature[i*increment+stnum-1]
		timeuse[i]=timing[i*increment+stnum-1]	
		killwaves tempwaveynorm,wcf,wsg
	endfor
	//timeuse=timeuse/60
	display peakuse vs timeuse
	appendtograph/R /C=(0,0,65535) tempuse vs timeuse
	Label left "Peak Wavelength (nm)";DelayUpdate
	Label bottom "Time, seconds";DelayUpdate
	Label right "Temperature, ¼C";DelayUpdate
	Legend/C/N=text0/F=0/H={30,4,10}/A=MC
	ModifyGraph grid(left)=1,grid(bottom)=1
	ModifyGraph mirror(bottom)=2
end

function trackgauss4a(wpre, stnum, increment, endnum, wlen,suffix, temperature,timing,start,stop,width) //with normalization, also plots peak vs temp
	string wpre,suffix
	wave temperature, timing,wlen
	variable stnum,endnum,increment, start, stop, width
	variable i,imax,maxloc,leftminloc,rightminloc,yleft,yright,slope,intcpt
	string tempwavey,peakname,errname,tempname,timename
	imax=floor((endnum-stnum)/increment)+1
	peakname="peak"+suffix
	errname="err"+suffix
	tempname="temp"+suffix
	timename="time"+suffix
	Make/N=(imax) /D $peakname
	Make/N=(imax) /D $errname
	Make/N=(imax) /D $tempname
	Make/N=(imax) /D $timename
	wave peakuse=$peakname
	wave erruse=$errname
	wave tempuse=$tempname
	wave timeuse=$timename
	for (i = 0; i < imax;i += 1)
		tempwavey=wpre+num2str(i*increment+stnum)
		wavestats/Q/R=[start,stop] $tempwavey
		maxloc=V_maxloc
		wavestats/Q/R=[maxloc-(1.5*width),maxloc-(0.5*width)] $tempwavey
		leftminloc=V_minloc
		yleft=V_min
		wavestats/Q/R=[maxloc+(0.5*width),maxloc+(1.5*width)] $tempwavey
		rightminloc=V_minloc
		yright=V_min
		slope=(yright-yleft)/(rightminloc-leftminloc)
		intcpt=yright-slope*rightminloc
		duplicate $tempwavey tempwaveynorm
		tempwaveynorm=tempwaveynorm/(slope*x+intcpt)
		CurveFit/Q/N=1/NTHR=0 gauss,  tempwaveynorm [maxloc-width/2, maxloc+width/2] /X=wlen /D
		wave wcf=W_coef
		wave wsg=W_sigma 
		peakuse[i]=wcf[2]
		erruse[i]=wsg[2]
		tempuse[i]=temperature[i*increment+stnum-1]
		timeuse[i]=timing[i*increment+stnum-1]	
		killwaves tempwaveynorm,wcf,wsg
	endfor
	//timeuse=timeuse/60
	display peakuse vs timeuse
	appendtograph/R /C=(0,0,65535) tempuse vs timeuse
	Label left "Peak Wavelength (nm)";DelayUpdate
	Label bottom "Time, seconds";DelayUpdate
	Label right "Temperature, degC";DelayUpdate
	Legend/C/N=text0/F=0/H={30,4,10}/A=MC
	ModifyGraph grid(left)=1,grid(bottom)=1
	ModifyGraph mirror(bottom)=2
	display peakuse vs tempuse
	Label left "Peak Wavelength (nm)";DelayUpdate
	Label bottom "Temperature, degC";DelayUpdate
	ModifyGraph grid(left)=1,grid(bottom)=1
	ModifyGraph mirror(bottom)=2
	ModifyGraph marker=8
end

function trackgauss5(wpre, stnum, increment, endnum, wlen,suffix, temperature,timing,start,stop,width) //with normalization after a poly[7] fit
	string wpre,suffix
	wave temperature, timing,wlen
	variable stnum,endnum,increment, start, stop, width
	variable i,imax,maxloc
	string tempwavey,peakname,errname,tempname,timename
	imax=floor((endnum-stnum)/increment)+1
	peakname="peak"+suffix
	errname="err"+suffix
	tempname="temp"+suffix
	timename="time"+suffix
	Make/N=(imax) /D $peakname
	Make/N=(imax) /D $errname
	Make/N=(imax) /D $tempname
	Make/N=(imax) /D $timename
	wave peakuse=$peakname
	wave erruse=$errname
	wave tempuse=$tempname
	wave timeuse=$timename
	wave W_fitConstants, W_coef
	for (i = 0; i < imax;i += 1)
		tempwavey=wpre+num2str(i*increment+stnum)
		CurveFit/Q/N=1/NTHR=0 poly_XOffset 7,  $tempwavey /X=wlen /D 
		wave normcoef=W_coef
		variable x0=W_fitConstants[0]
		duplicate/O $tempwavey tempwaveynorm
		tempwaveynorm=tempwaveynorm/(normcoef[0]+normcoef[1]*(wlen-x0)+normcoef[2]*(wlen-x0)^2+normcoef[3]*(wlen-x0)^3+normcoef[4]*(wlen-x0)^4+normcoef[5]*(wlen-x0)^5+normcoef[6]*(wlen-x0)^6)
		wavestats/Q/R=[start,stop] tempwaveynorm
		maxloc=V_maxloc
		CurveFit/Q/N=1/NTHR=0 gauss,  tempwaveynorm [maxloc-width/2, maxloc+width/2] /X=wlen /D
		wave wcf=W_coef
		wave wsg=W_sigma 
		peakuse[i]=wcf[2]
		erruse[i]=wsg[2]
		tempuse[i]=temperature[i*increment+stnum-1]
		timeuse[i]=timing[i*increment+stnum-1]	
		killwaves tempwaveynorm,wcf,wsg
	endfor
	//timeuse=timeuse/60
	display peakuse vs timeuse
	appendtograph/R /C=(0,0,65535) tempuse vs timeuse
	Label left "Peak Wavelength (nm)";DelayUpdate
	Label bottom "Time, seconds";DelayUpdate
	Label right "Temperature, ¼C";DelayUpdate
	Legend/C/N=text0/F=0/H={30,4,10}/A=MC
	ModifyGraph grid(left)=1,grid(bottom)=1
	ModifyGraph mirror(bottom)=2
end

function trackgaussall(wpre,stnum,increment,endnum, wlen,suffix,temperature, timing,peak1,period,nopeaks,width) //analyzes all peaks
	string wpre,suffix
	wave temperature, timing,wlen
	variable peak1, period,nopeaks,stnum,endnum,increment, width
	variable i
	string suffixi
	for  (i = 0; i < nopeaks;i += 1)
	suffixi=suffix+num2str(i)
	trackgauss4(wpre, stnum, increment, endnum, wlen,suffixi, temperature,timing,peak1+i*period-0.5*width,peak1+i*period+0.5*width,width)
	endfor
end

function trackgaussall2(wpre,stnum,increment,endnum, wlen,suffix,temperature, timing,peakswave,width,shiftmax) //analyzes all peaks
	string wpre,suffix
	wave temperature, timing,wlen,peakswave
	variable stnum,endnum,increment, width,shiftmax
	variable i,len
	string suffixi
	wavestats/Q peakswave
	len= V_npnts
	for  (i = 0; i < len;i += 1)
	suffixi=suffix+num2str(i)
	trackgauss4(wpre, stnum, increment, endnum, wlen,suffixi, temperature,timing,peakswave[i]-round(0.2*width),peakswave[i]+round(0.2*width+shiftmax),width)
	endfor
end

function trackgaussall3(wpre,stnum,increment,endnum, wlen,suffix,temperature, timing,peaksinitial,peaksfinal,width) //analyzes all peaks
	string wpre,suffix
	wave temperature, timing,wlen,peaksinitial,peaksfinal
	variable stnum,endnum,increment, width
	variable i,len
	string suffixi
	wavestats/Q peaksinitial
	len= V_npnts
	for  (i = 0; i < len;i += 1)
	suffixi=suffix+num2str(i)
	trackgauss4(wpre, stnum, increment, endnum, wlen,suffixi, temperature,timing,peaksinitial[i]-round(0.2*width),peaksfinal[i]+round(0.2*width),width)
	endfor
end

function peakavg(peakwave,riwave,stravg,stpavg, suffix)
	wave peakwave,riwave,stravg,stpavg
	string suffix
	
	variable i,len
	string avguse,erruse
	avguse="peakavg"+suffix
	erruse="peakerr"+suffix
	duplicate riwave $avguse,$erruse
	wave avgusem=$avguse
	wave errusem=$erruse
	wavestats/Q riwave
	len= V_npnts
	for  (i = 0; i < len;i += 1)
	wavestats/Q/R=[stravg[i],stpavg[i]] peakwave
	avgusem[i]=V_avg
	errusem[i]=V_sdev
	endfor
	display avgusem vs riwave
	ModifyGraph mode=3
	//ErrorBars avgusem Y,wave=(errusem,errusem)
	CurveFit/NTHR=0 line  avgusem /X=riwave /D
	//addslope to graph
end


function allpeakavg(nopeaks,riwave,stravg,stpavg, suffix)
	string suffix
	variable nopeaks
	wave riwave,stravg,stpavg
	
	variable i
	string peakw,sensoutstr,sensouterrstr,wlenoutstr,suffixi
	wave W_coef, W_sigma
	sensoutstr="sensout"+suffix
	sensouterrstr="sensouterr"+suffix
	wlenoutstr="wlenout"+suffix
	Make/N=(nopeaks) /D $sensoutstr
	Make/N=(nopeaks) /D $sensouterrstr
	Make/N=(nopeaks) /D $wlenoutstr
	wave sensout=$sensoutstr
	wave sensouterr=$sensouterrstr
	wave wlenout=$wlenoutstr
	
	for  (i = 0; i < nopeaks;i += 1)
	peakw="peak"+suffix+num2str(i)
	suffixi=suffix+num2str(i)
	wave peakwave=$peakw
	peakavg(peakwave,riwave,stravg,stpavg, suffixi)
	sensout[i]=W_coef[1]
	sensouterr[i]=W_sigma[1]
	wavestats/Q peakwave
	wlenout[i]=V_avg
	killwaves peakwave
	endfor
	display sensout vs wlenout
end

function peakstats(nopeaks,stravg,stpavg, suffix)
	string suffix
	variable nopeaks,stravg,stpavg
	
	variable i
	string peakw,wlenoutstr,suffixi
	wave W_coef, W_sigma
	wlenoutstr="wlenout"+suffix
	Make/N=(nopeaks) /D $wlenoutstr
	wave wlenout=$wlenoutstr
	for  (i = 0; i < nopeaks;i += 1)
	peakw="peak"+suffix+num2str(i)
	suffixi=suffix+num2str(i)
	wave peakwave=$peakw 
	wavestats/Q/R=[stravg,stpavg] peakwave
	wlenout[i]=V_avg
	killwaves peakwave
	endfor
	display sensout vs wlenout
end

function allpeakfits(wpre,nopeaks,riwave,suffix)
	string wpre,suffix
	variable nopeaks
	wave riwave
	
	variable i
	string peakw,sensoutstr,sensouterrstr,wlenoutstr
	wave W_coef, W_sigma
	sensoutstr="sensout"+suffix
	sensouterrstr="sensouterr"+suffix
	wlenoutstr="wlenout"+suffix
	Make/N=(nopeaks) /D $sensoutstr
	Make/N=(nopeaks) /D $sensouterrstr
	Make/N=(nopeaks) /D $wlenoutstr
	wave sensout=$sensoutstr
	wave sensouterr=$sensouterrstr
	wave wlenout=$wlenoutstr
	
	for  (i = 0; i < nopeaks;i += 1)
	peakw=wpre+num2str(i)
	wave peakwave=$peakw
	CurveFit/Q/NTHR=0 line  peakwave /X=riwave /D
	sensout[i]=W_coef[1]
	sensouterr[i]=W_sigma[1]
	wavestats/Q peakwave
	wlenout[i]=V_avg
	killwaves peakwave
	endfor
	display sensout vs wlenout
end

function justgaussfit(wpre, stnum, increment, endnum, wlen,suffix, temperature,timing,start,stop,width) //with normalization
	string wpre,suffix
	wave temperature, timing,wlen
	variable stnum,endnum,increment,width,start,stop
	variable i,imax,maxloc,leftminloc,rightminloc,yleft,yright,slope,intcpt
	string tempwavey,peakname,errname,tempname,timename
	imax=floor((endnum-stnum)/increment)+1
	peakname="peak"+suffix
	errname="err"+suffix
	tempname="temp"+suffix
	timename="time"+suffix
	Make/N=(imax) /D $peakname
	Make/N=(imax) /D $errname
	Make/N=(imax) /D $tempname
	Make/N=(imax) /D $timename
	wave peakuse=$peakname
	wave erruse=$errname
	wave tempuse=$tempname
	wave timeuse=$timename
	for (i = 0; i < imax;i += 1)
		tempwavey=wpre+num2str(i*increment+stnum)
		wavestats/Q/R=[start,stop] $tempwavey
		maxloc=V_maxloc
		wavestats/Q/R=[maxloc-(1.5*width),maxloc-(0.5*width)] $tempwavey
		leftminloc=V_minloc
		yleft=V_min
		wavestats/Q/R=[maxloc+(0.5*width),maxloc+(1.5*width)] $tempwavey
		rightminloc=V_minloc
		yright=V_min
		slope=(yright-yleft)/(rightminloc-leftminloc)
		intcpt=yright-slope*rightminloc
		duplicate $tempwavey tempwaveynorm
		tempwaveynorm=tempwaveynorm/(slope*x+intcpt)
		wavestats/Q/R=[start,stop] tempwaveynorm
		maxloc=V_maxloc
		CurveFit/Q/N=1/NTHR=0 gauss,  tempwaveynorm [maxloc-width/2, maxloc+width/2] /X=wlen /D
		wave wcf=W_coef
		wave wsg=W_sigma 
		peakuse[i]=wcf[2]
		erruse[i]=wsg[2]
		tempuse[i]=temperature[i*increment+stnum-1]
		timeuse[i]=timing[i*increment+stnum-1]	
		killwaves tempwaveynorm,wcf,wsg
	endfor
	//timeuse=timeuse/60
	display peakuse vs timeuse
	appendtograph/R /C=(0,0,65535) tempuse vs timeuse
	Label left "Peak Wavelength (nm)";DelayUpdate
	Label bottom "Time, seconds";DelayUpdate
	Label right "Temperature, ¼C";DelayUpdate
	Legend/C/N=text0/F=0/H={30,4,10}/A=MC
	ModifyGraph grid(left)=1,grid(bottom)=1
	ModifyGraph mirror(bottom)=2
end
