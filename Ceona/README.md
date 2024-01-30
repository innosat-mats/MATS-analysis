Aurora Analysis of the MATS satellite

# Insert link to report

--- Project Description ---
The scripts in this folder are a part of a master thesis work. They can be used to make a deeper aurora analysis of aurora observed in the images taken by the MATS satellite.


--- Table of Contents ---

1. Keogram.py
Contains functions to create centerstrip objects and corresponding matrix (Keogram) from a list of ccd objects.
    - class CenterStrip
    - def makeStripMatrix()
    - def getTPLatitudes(CCDs)
    - def getSatDates(CCDs)
    - def get_stripRow(CCDs)
2. aurora_analysis.py
Contains functions to calculate aurora parameters, such as location and intensity.
    - def TPpos(ccd)
    - def save_TPMLT(CCDs, filename)
    - def set_strip_spec(strip,ccd)
    - def col_pos(ccd, column)
    - def set_aurora_spec(strip,ccd,row)
    - def IntensityPeak(strip)
    - def get_all_altitudes(strips)
    - def get_aurora_max(aurorastrips,filedate)
    - def save_strips(strips,filename,structname)
3. Auroradetecion.py
Functions to create and save list of strip objects from a list of CCDs. And a separate list of strips with aurora observed. This detection algorithm uses polynomial regression to remove a background gradient for easer detection.
    - def gradientmatrix(ccds..)
    - def get_strips(items,numdays,filedate)
4. gradientOverview.py
Create pdf overviews with Keogram plots for every orbit and hemisphere, the code creates keograms with the background gradient reduced from polynomial regression. Added red dots in the plots for the aurora strips maximum intensity points. 
5. orbitoverview.py
Create pdf overviews with Keogram plots for every orbit and hemisphere. 

--- How to Install and Run the Project ---
Python 3.9 or 3.10
Skyfield
Python wrapper AACGMV2

Other needed repositories are: MATS-L1-processing and MATS-utility-functions, can be found on github "innosat-Mats".

--- How to Use the Project ---
Use the Runfile.py to read in raw MATS data with chosen filters, preferably to save as a pickle, and to run the wanted functions from the other scripts.

--- Include Credits ---
MISU - Stockholm University
Nickolay Ivchenko - Assistant professor KTH


--- License ---
MIT License

Copyright (c) 2023 innosat-mats

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.


