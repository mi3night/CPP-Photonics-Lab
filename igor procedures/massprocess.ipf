#pragma rtGlobals=1		// Use modern global access method.
function massplot(wpre, stnum, increment, endnum, wlen)
	string wpre
	wave wlen
	variable stnum,endnum,increment
	variable i,imax
	string tempwave,tempwave1
	imax=floor((endnum-stnum)/increment)+1
	tempwave1=wpre+num2str(stnum)
	wave tempwave2=$tempwave1
	display tempwave2 vs wlen
	for (i = 0; i < imax;i += 1)
		tempwave=wpre+num2str(i*increment+stnum)
		appendtograph $tempwave vs wlen
	endfor
	Label left "Optical Power, mW";DelayUpdate
	Label bottom "Wavelength, nm";DelayUpdate
	ModifyGraph grid(left)=1,grid(bottom)=1
	ModifyGraph mirror=2
end

function masskill(winpre, stnum, increment, endnum) //kills all the windows instantly
	string winpre
	variable stnum,endnum,increment
	variable i,imax
	string tempname
	imax=floor((endnum-stnum)/increment)+1
	for (i = 0; i < imax;i += 1)	
		tempname=winpre+num2str(i*increment+stnum)
		KillWindow $tempname
	endfor
end

function masskillwaves(wpre, stnum, increment, endnum) //kills all the intended waves instantly
	string wpre
	variable stnum,endnum,increment
	variable i,imax
	string tempname
	imax=floor((endnum-stnum)/increment)+1
	for (i = 0; i < imax;i += 1)	
		tempname=wpre+num2str(i*increment+stnum)
		KillWaves $tempname
	endfor
end


function moviemaker(moviename,wpre, stnum, increment, endnum, wlen,timing,temp,ymin,ymax) //makes a movie from data
	string wpre,moviename
	wave wlen,timing,temp
	variable stnum,endnum,increment,ymin,ymax
	variable i,imax
	string tempwave,annotate1
	imax=floor((endnum-stnum)/increment)+1
	tempwave=wpre+num2str(stnum)
	wave tempwave1=$tempwave
	duplicate tempwave1 runwave
	display/K=1 /N=tempgr runwave vs wlen
	Label left "Optical Power, mW";DelayUpdate
	Label bottom "Wavelength, nm";DelayUpdate
	SetAxis left ymin,ymax
	ModifyGraph grid(left)=1,grid(bottom)=1
	ModifyGraph mirror=2
	NewMovie as moviename
	for (i = 0; i < imax;i += 1)
		tempwave=wpre+num2str(i*increment+stnum)
		wave tempwave1=$tempwave
		runwave=tempwave1
		annotate1="time="+num2str(round(timing[i*increment+stnum-1]))+" sec, T="+num2str(temp[i*increment+stnum-1])+"\So\MC"
		TextBox/W=tempgr /A=LT /F=0 annotate1
		DoUpdate
		AddMovieFrame
	endfor
	CloseMovie
	KillWindow tempgr
	Killwaves runwave
end