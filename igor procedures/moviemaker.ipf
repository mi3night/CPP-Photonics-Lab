#pragma rtGlobals=1		// Use modern global access method.
function moviemaker(moviename,wpre, stnum, increment, endnum, wlen,timing,temper,xmin,xmax,ymin,ymax)
	string wpre,moviename
	wave wlen,timing,temper
	variable stnum,endnum,increment,xmin,xmax,ymin,ymax
	variable i,imax
	string tempwave,annotate1,annotate2
	imax=floor((endnum-stnum)/increment)+1
	tempwave=wpre+num2str(stnum)
	wave tempwave1=$tempwave
	duplicate tempwave1 runwave
	display/K=1 /N=tempgr runwave vs wlen
	appendtograph tempwave1 vs wlen
	ModifyGraph rgb($tempwave)=(0,0,0)
	Label left "Optical Power, mW";DelayUpdate
	Label bottom "Wavelength, nm";DelayUpdate
	SetAxis bottom xmin,xmax
	SetAxis left ymin,ymax
	ModifyGraph grid(left)=1,grid(bottom)=1
	ModifyGraph mirror=2
	NewMovie as moviename
	for (i = 0; i < imax;i += 1)
		tempwave=wpre+num2str(i*increment+stnum)
		wave tempwave1=$tempwave
		runwave=tempwave1
		annotate1="time="+num2str(round(timing[i*increment+stnum-1]))+" sec"
		annotate2="T="+num2str(round(temper[i*increment+stnum-1]*100)/100)+" degC"
		TextBox/W=tempgr /A=LT /F=0 annotate1
		TextBox/W=tempgr /A=RT /F=0 annotate2
		DoUpdate
		AddMovieFrame
	endfor
	CloseMovie
	KillWindow tempgr
	Killwaves runwave
end