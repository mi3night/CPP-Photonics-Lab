#pragma rtGlobals=1		// Use modern global access method.

function normalize(wavein,wlen)
	wave wavein,wlen
	CurveFit/Q/N=1/NTHR=0 poly_XOffset 7,  wavein /X=wlen /D 
	wave normcoef=W_coef
	wave W_fitConstants
	variable x0=W_fitConstants[0]
	return wavein/(normcoef[0]+normcoef[1]*(wlen-x0)+normcoef[2]*(wlen-x0)^2+normcoef[3]*(wlen-x0)^3+normcoef[4]*(wlen-x0)^4+normcoef[5]*(wlen-x0)^5+normcoef[6]*(wlen-x0)^6)
end

function normalize1(wavein,wlen,coef) // coef is W_coef after poly_offset fit to the first wave in the sequence, 7th entry is the constant for poly_offset
	wave wavein,wlen,coef
	variable x0
	x0=coef[7]
	return wavein/(coef[0]+coef[1]*(wlen-x0)+coef[2]*(wlen-x0)^2+coef[3]*(wlen-x0)^3+coef[4]*(wlen-x0)^4+coef[5]*(wlen-x0)^5+coef[6]*(wlen-x0)^6)
end

function normalizeref(wavein,ref) //normalizes to a reference wave
	wave wavein,ref
	return wavein/ref
end

function normalizeall(wpre, stnum, increment, endnum, wlen)
	string wpre
	wave wlen
	variable stnum, endnum, increment
	string tempwavey, nwave
	variable i, imax
	imax=floor((endnum-stnum)/increment)+1
	for (i = 0; i < imax;i += 1)
		tempwavey=wpre+num2str(i*increment+stnum)
		nwave="n"+wpre+num2str(i*increment+stnum)
		Make/N=1001 /D $nwave
		wave normwave=$nwave
		normwave=normalize($tempwavey,wlen)
	endfor
end

function normalizeall1(npre,wpre, stnum, increment, endnum, wlen,coef)
	string wpre,npre
	wave wlen,coef
	variable stnum, endnum, increment
	string tempwavey, nwave
	variable i, imax
	imax=floor((endnum-stnum)/increment)+1
	for (i = 0; i < imax;i += 1)
		tempwavey=wpre+num2str(i*increment+stnum)
		nwave=npre+wpre+num2str(i*increment+stnum)
		Make/N=1001 /D $nwave
		wave normwave=$nwave
		normwave=normalize1($tempwavey,wlen,coef)
	endfor
end

function normalizeallref(npre,wpre, stnum, increment, endnum, ref)
	string wpre,npre
	wave ref
	variable stnum, endnum, increment
	string tempwavey, nwave
	variable i, imax
	imax=floor((endnum-stnum)/increment)+1
	for (i = 0; i < imax;i += 1)
		tempwavey=wpre+num2str(i*increment+stnum)
		nwave=npre+num2str(i*increment+stnum)
		Make/N=1001 /D $nwave
		wave normwave=$nwave
		normwave=normalizeref($tempwavey,ref)
	endfor
end