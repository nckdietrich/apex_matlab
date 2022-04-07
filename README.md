# apex_matlab
Matlab wrapping code for apex implemented in Fortran. Check matlab functions for inputs and outputs. 


To compile apex:
    ifort -o apex apex.f90 test.f

(Note for myself)
To compile apex calling fortran code:
 - Compile all files individually to make .o outputs
    ifort -c apex.f90 matlab_apex_mall.f90
 - Compile the executable
    ifort apex.o matlab_apex_mall.o -o matlab_apex_mall
 - Run executable
    ./matlab_apex_mall

Make sure executable functions are in the top level folder

Available Apex commands:
 - apex_mall: Convert from geodetic coordinates to apex/quasi-dipole coordinates + magnetic fields + basis vectors
 - apex_q2g: Convert from quasi-dipole to geodetic coordinates
 - apex_m2g: Convert from modified apex coordinates to geodetic coordinates
