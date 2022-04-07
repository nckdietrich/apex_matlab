# apex_matlab
Matlab wrapping code for apex implemented in Fortran. Check matlab functions for inputs and outputs. 


To compile apex:
    ifort -o apex apex.f90 test.f


To compile apex calling fortran code:
 - Compile all files individually to make .o outputs
    ifort -c apex.f90 matlab_apex_mall.f90
 - Compile the executable
    ifort apex.o matlab_apex_mall.o -o a.out
 - Run executable
    ./a.out

Make sure executable functions are in the top level folder


Apex commands:
 - apex_mall: Convert from geodetic coordinates to apex/quasi-dipole coordinates
 - apex_q2g: Convert from quasi-dipole to geodetic coordinates
 - apex_m2g: Convert from modified apex coordinates to geodetic coordinates
