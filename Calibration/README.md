## PMT equalization and gain calibration

Usage Instructions:
1. root -l 'drawIntADC_3x8.C(run#, # of events, -1, -1, "save")' // for saving

    or root -l 'drawIntADC_3x8.C(run#, # of events, 1, 1, "draw"")' // for drawing

2. root -l 'Fit_HV.C(row#)' 

File Descriptions:

drawIntADC_3x8.C: Draws or saves IntADC (integrated ADC) plots.

Fit_HV.C: Calculates the appropriate high voltage (HV) using fitting.
