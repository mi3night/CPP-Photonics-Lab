#pragma rtGlobals=3		// Use modern global access method and strict wave access.
function findmax(w,maxpeaks,start,stop)
Wave w  // peak data
Variable maxPeaks,start,stop
Make/O/N=(maxPeaks) peakPositionsX= NaN, peakPositionsY= NaN    
Variable peaksFound=0
Variable startP=start
Variable endP= stop
do
    FindPeak/B=2/I/P/Q/R=[startP,endP] w
    // FindPeak outputs are V_Flag, V_PeakLoc, V_LeadingEdgeLoc,
    // V_TrailingEdgeLoc, V_PeakVal, and V_PeakWidth.
    if( V_Flag != 0 )
        break
    endif
   
    peakPositionsX[peaksFound]=pnt2x(w,V_PeakLoc)
    peakPositionsY[peaksFound]=V_PeakVal
    peaksFound += 1
   
    startP= V_TrailingEdgeLoc+1
while( peaksFound < maxPeaks )

if( peaksFound )
    Redimension/N=(peaksFound) peakPositionsX, peakPositionsY

    DoWindow/F ShowPeaks
    if(V_Flag == 0 )
        Display/N=ShowPeaks w
        AppendToGraph/W=ShowPeaks peakPositionsY vs peakPositionsX
        ModifyGraph/W=ShowPeaks mode(peakPositionsY)=3,marker(peakPositionsY)=19, rgb(peakPositionsY)=(0,0,65535)
    endif
else
    DoAlert 0, "No peaks found"
    KillWaves/Z peakPositionsX, peakPositionsY
endif
Edit peakPositionsX, peakPositionsY
end

function findmin(w,maxpeaks,start,stop)
Wave w  // peak data
Variable maxPeaks,start,stop
Make/O/N=(maxPeaks) peakPositionsX= NaN, peakPositionsY= NaN    
Variable peaksFound=0
Variable startP=start
Variable endP= stop
do
    FindPeak/N/B=2/I/P/Q/R=[startP,endP] w
    // FindPeak outputs are V_Flag, V_PeakLoc, V_LeadingEdgeLoc,
    // V_TrailingEdgeLoc, V_PeakVal, and V_PeakWidth.
    if( V_Flag != 0 )
        break
    endif
   
    peakPositionsX[peaksFound]=pnt2x(w,V_PeakLoc)
    peakPositionsY[peaksFound]=V_PeakVal
    peaksFound += 1
   
    startP= V_TrailingEdgeLoc+1
while( peaksFound < maxPeaks )

if( peaksFound )
    Redimension/N=(peaksFound) peakPositionsX, peakPositionsY

    DoWindow/F ShowPeaks
    if(V_Flag == 0 )
        Display/N=ShowPeaks w
        AppendToGraph/W=ShowPeaks peakPositionsY vs peakPositionsX
        ModifyGraph/W=ShowPeaks mode(peakPositionsY)=3,marker(peakPositionsY)=19, rgb(peakPositionsY)=(0,0,65535)
    endif
else
    DoAlert 0, "No peaks found"
    KillWaves/Z peakPositionsX, peakPositionsY
endif
Edit peakPositionsX, peakPositionsY
end