# LPRM
LPRM: user-friendly iteration-free combined Local Polynomial and Rational method for measurements with multiple inputs

   LPRM: iteration free Local Polynomial/Rational Method
   Copyright (C) 2010  Dr. Péter Zoltán Csurcsia
   Details on license can be found in license.txt

   Start by typing 'help LPRM'

   LPRM.m contains the FRM generation methods 

   demo.m contains a serios of demos that is recommended to start with

   test_data.mat contains a measurement of a car frame testing using multisine generated
   by the MUMI toolbox
	  
   Sources folder contains:
    - the default settings file (Test_settings.m) that can be modified
   and auxiliary codes:
   - Estimate_N.m: estimates the length of the period    
   - Estimate_N.m: estimates the number of periods
   - Estimate_Ptr.m: estimates the number of transient (delay) blocks
   - checkCorrelation.m: checks the correlation between (input) channels
   - checkDefaultValuesAIO: initialize the options structure and estimate parameters
   - segmentData.m: segments data in the time and frequency domain
   - generateDecadeBand.m and generateOctaveBand.m that generate non-decade frequency lines
   - determineSolver.m: determines which solver has to be used
   - lrm_fd.m: it is the LPM/LRM code


   Related article with technical details on the recommended nonlinear framework can be found in
   manuscript.pdf, for scientific reference, please use:
        Péter Zoltán Csurcsia, Bart Peeters, Johan Schoukens:
        User-friendly nonlinear nonparametric estimation framework for vibro-acoustic industrial measurements with multiple inputs,
        Mechanical Systems and Signal Processing, Volume 145, 2020
        https://doi.org/10.1016/j.ymssp.2020.106926.
