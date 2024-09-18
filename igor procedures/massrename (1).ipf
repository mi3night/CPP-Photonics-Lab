#pragma rtGlobals=1		// Use modern global access method.

function massrename(wprein, stnum,endnum, wpreout)
	string wprein, wpreout
	variable stnum,endnum
	variable i,imax
	string wavein,waveout
	imax=endnum-stnum+1
	for (i = 0; i < imax;i += 1)
		wavein=wprein+num2str(i+stnum)
		waveout=wpreout+num2str(i)
		rename $wavein $waveout
	endfor
end

function massrename2(wprein, stnum,endnum, wpreout,startout)
	string wprein, wpreout
	variable stnum,endnum,startout
	variable i,imax
	string wavein,waveout
	imax=endnum-stnum+1
	for (i = 0; i < imax;i += 1)
		wavein=wprein+num2str(i+stnum)
		waveout=wpreout+num2str(i+startout)
		rename $wavein $waveout
	endfor
end