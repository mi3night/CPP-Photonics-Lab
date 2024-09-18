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

function massappend(wpre, stnum, increment, endnum, wlen)
	string wpre
	wave wlen
	variable stnum,endnum,increment
	variable i,imax
	string tempwave,tempwave1
	imax=floor((endnum-stnum)/increment)+1
	tempwave1=wpre+num2str(stnum)
	for (i = 0; i < imax;i += 1)
		tempwave=wpre+num2str(i*increment+stnum)
		appendtograph $tempwave vs wlen
	endfor
end

function masskill(winpre, stnum, increment, endnum)
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
