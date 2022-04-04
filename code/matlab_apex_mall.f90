! Running apex_mall - Converts geodetic coords to apex/quasi-dipole coords
! Author: Nick Dietrich
! Version 4/1/22

program matlab_apex_mall

use apex
use netcdf

implicit none

!! Declaring 
character(len=20) :: file_in = 'geodetic_pos.nc'
character(len=128) :: file_out = 'apex_pos.nc'
integer :: in_ncid, dimid, varid
integer :: out_ncid, dimid_nloc, dimid_locdim, dimid_basisdim, dimid_enudim
integer :: nloc, loc_ind, ier
real :: date, altmax, hr
real, allocatable :: loc(:,:), glons(:), glats(:), alts(:)
! Apex output variables
real :: bmag,si,alon,malat,vmp,W,D,Be3,sim,qdlat,F
real :: B(3),bhat(3),d1(3),d2(3),d3(3),e1(3),e2(3),e3(3),f1(3),f2(3) &
 ,f3(3),g1(3),g2(3),g3(3)
real, allocatable :: apex_loc(:,:), qdipole_loc(:,:)
real, allocatable :: mag_field_comp(:,:), mag_field_vector(:,:), bmagitude(:), &
                        sine_mag_inc(:), mag_potential(:), W_richmond(:), &
                        D_richmond(:), be3_richmond(:), sinIm_richmond(:), &
                        F_richmond(:), ier_list(:)
real, allocatable :: d_basis(:,:,:), e_basis(:,:,:), f_basis(:,:,:), g_basis(:,:,:)

!! Running
print *, 'Converting Geodetic Coords to Apex/Quasi-Dipole Coords'
! Read in information
! Get vector length
call check( nf90_open(file_in, NF90_NOWRITE, in_ncid) )
call check( nf90_inq_dimid(in_ncid, 'nloc', dimid))
call check( nf90_inquire_dimension(in_ncid, dimid, len=nloc))
allocate(loc(nloc,3))
! Locations [lon, lat, alt]
call check( nf90_inq_varid(in_ncid, 'loc', varid) )
call check( nf90_get_var(in_ncid, varid, values=loc) )
! date
call check( nf90_inq_varid(in_ncid, 'date', varid) )
call check( nf90_get_var(in_ncid, varid ,values=date))
! close
call check( nf90_close(in_ncid) )

! Running apex code
hr = 90 ! Typical value used in TIEGCM [km]
altmax = MAXVAL(loc(:,3))
allocate(glons(nloc))
allocate(glats(nloc))
allocate(alts(nloc))
glons = loc(:,1)
glats = loc(:,2)
alts = loc(:,3)

! Set up interp arrays (takes longest to run) + allocate
call apex_setup(date, altmax)
allocate(apex_loc(nloc,3))
allocate(qdipole_loc(nloc,3))
allocate(mag_field_comp(nloc,3))
allocate(mag_field_vector(nloc,3))
allocate(bmagitude(nloc))
allocate(sine_mag_inc(nloc))
allocate(mag_potential(nloc))
allocate(W_richmond(nloc))
allocate(D_richmond(nloc))
allocate(be3_richmond(nloc))
allocate(sinIm_richmond(nloc))
allocate(F_richmond(nloc))
allocate(d_basis(nloc,3,3))
allocate(e_basis(nloc,3,3))
allocate(f_basis(nloc,3,3))
allocate(g_basis(nloc,3,3))
allocate(ier_list(nloc))
! Run apex
print *, 'Looping through locations...'
do loc_ind = 1, nloc
    call apex_mall( glats(loc_ind),glons(loc_ind),alts(loc_ind),hr, &
        B,bhat,bmag,si,alon,malat,vmp,W,D,Be3,sim,d1,d2,d3,e1,e2,e3,&
        qdlat,F,f1,f2,f3,g1,g2,g3,ier )

    ! Save to variables
    apex_loc(loc_ind,:) = [alon, malat, alts(loc_ind)]
    qdipole_loc(loc_ind,:) = [alon, qdlat, alts(loc_ind)]
    mag_field_comp(loc_ind,:) = B
    mag_field_vector(loc_ind,:) = bhat
    bmagitude(loc_ind) = bmag
    sine_mag_inc(loc_ind) = si
    mag_potential(loc_ind) = vmp
    W_richmond(loc_ind) = W
    D_richmond(loc_ind) = D
    be3_richmond(loc_ind) = Be3
    sinIm_richmond(loc_ind) = sim
    F_richmond(loc_ind) = F
    d_basis(loc_ind,1,:) = d1
    d_basis(loc_ind,2,:) = d2
    d_basis(loc_ind,3,:) = d3
    e_basis(loc_ind,1,:) = e1
    e_basis(loc_ind,2,:) = e2
    e_basis(loc_ind,3,:) = e3
    f_basis(loc_ind,1,:) = f1
    f_basis(loc_ind,2,:) = f2
    f_basis(loc_ind,3,:) = f3
    g_basis(loc_ind,1,:) = g1
    g_basis(loc_ind,2,:) = g2
    g_basis(loc_ind,3,:) = g3
    ier_list(loc_ind) = ier
enddo

! Write out data (saving everyting)
print *, 'Writing out to file...'
call check( nf90_create(file_out, NF90_NETCDF4, out_ncid) )
! Define dimensions
call check( nf90_def_dim(out_ncid, 'nloc', nloc, dimid_nloc) )
call check( nf90_def_dim(out_ncid, 'locDim', 3, dimid_locdim))
call check( nf90_def_dim(out_ncid, 'basisDim', 3, dimid_basisdim))
call check( nf90_def_dim(out_ncid, 'ENUDim', 3, dimid_enudim))

! Writing out variables
call check( nf90_def_var(out_ncid, 'apex_loc', NF90_DOUBLE, [dimid_nloc, dimid_locdim], varid) )
call check( nf90_put_var(out_ncid, varid, apex_loc) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'Modified Apex Coordinates') )
call check( nf90_put_att(out_ncid, varid, 'arrangement', '[lon, lat, alt]') )
call check( nf90_put_att(out_ncid, varid, 'units', '[deg, deg, km]') )

call check( nf90_def_var(out_ncid, 'qdipole_loc', NF90_DOUBLE, [dimid_nloc, dimid_locdim], varid) )
call check( nf90_put_var(out_ncid, varid, qdipole_loc) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'Quasi-dipole Coordinates') )
call check( nf90_put_att(out_ncid, varid, 'arrangement', '[lon, lat, alt]') )
call check( nf90_put_att(out_ncid, varid, 'units', '[deg, deg, km]') )

call check( nf90_def_var(out_ncid, 'mag_field_comp', NF90_DOUBLE, [dimid_nloc, dimid_enudim], varid) )
call check( nf90_put_var(out_ncid, varid, mag_field_comp) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'Magnetic field components') )
call check( nf90_put_att(out_ncid, varid, 'arrangement', '[east, north, up]') )
call check( nf90_put_att(out_ncid, varid, 'units', '[nT]') )

call check( nf90_def_var(out_ncid, 'mag_field_vector', NF90_DOUBLE, [dimid_nloc, dimid_enudim], varid) )
call check( nf90_put_var(out_ncid, varid, mag_field_vector) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'Components of unit vector along geomagnetic field direction') )
call check( nf90_put_att(out_ncid, varid, 'arrangement', '[east, north, up]') )

call check( nf90_def_var(out_ncid, 'bmagitude', NF90_DOUBLE, [dimid_nloc], varid) )
call check( nf90_put_var(out_ncid, varid, bmagitude) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'Magnitude of magnetic field') )
call check( nf90_put_att(out_ncid, varid, 'units', '[nT]') )

call check( nf90_def_var(out_ncid, 'sine_mag_inc', NF90_DOUBLE, [dimid_nloc], varid) )
call check( nf90_put_var(out_ncid, varid, sine_mag_inc) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'Sine of magnetic inclination') )

call check( nf90_def_var(out_ncid, 'mag_potential', NF90_DOUBLE, [dimid_nloc], varid) )
call check( nf90_put_var(out_ncid, varid, mag_potential) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'Magnetic potential') )
call check( nf90_put_att(out_ncid, varid, 'units', '[T.m]') )

call check( nf90_def_var(out_ncid, 'W_richmond', NF90_DOUBLE, [dimid_nloc], varid) )
call check( nf90_put_var(out_ncid, varid, W_richmond) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'W of Richmond [1995] reference') )
call check( nf90_put_att(out_ncid, varid, 'units', '[km^2/nT]') )

call check( nf90_def_var(out_ncid, 'D_richmond', NF90_DOUBLE, [dimid_nloc], varid) )
call check( nf90_put_var(out_ncid, varid, D_richmond) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'D of Richmond [1995]') )

call check( nf90_def_var(out_ncid, 'be3_richmond', NF90_DOUBLE, [dimid_nloc], varid) )
call check( nf90_put_var(out_ncid, varid, be3_richmond) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'B_e3 of Richmond [1995] (= Bmag/D)') )
call check( nf90_put_att(out_ncid, varid, 'units', '[nT]') )

call check( nf90_def_var(out_ncid, 'sinIm_richmond', NF90_DOUBLE, [dimid_nloc], varid) )
call check( nf90_put_var(out_ncid, varid, sinIm_richmond) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'sin(I_m) described in Richmond [1995]') )

call check( nf90_def_var(out_ncid, 'F_richmond', NF90_DOUBLE, [dimid_nloc], varid) )
call check( nf90_put_var(out_ncid, varid, F_richmond) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'F described in Richmond [1995] for quasi-dipole coords') )

call check( nf90_def_var(out_ncid, 'd_basis', NF90_DOUBLE, [dimid_nloc, dimid_basisdim, dimid_enudim], varid) )
call check( nf90_put_var(out_ncid, varid, d_basis) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'd basis vectors') )
call check( nf90_put_att(out_ncid, varid, 'basis arrangement', '[east, north, up]') )

call check( nf90_def_var(out_ncid, 'e_basis', NF90_DOUBLE, [dimid_nloc, dimid_basisdim, dimid_enudim], varid) )
call check( nf90_put_var(out_ncid, varid, e_basis) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'e basis vectors') )
call check( nf90_put_att(out_ncid, varid, 'basis arrangement', '[east, north, up]') )

call check( nf90_def_var(out_ncid, 'f_basis', NF90_DOUBLE, [dimid_nloc, dimid_basisdim, dimid_enudim], varid) )
call check( nf90_put_var(out_ncid, varid, f_basis) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'f basis vectors') )
call check( nf90_put_att(out_ncid, varid, 'basis arrangement', '[east, north, up]') )

call check( nf90_def_var(out_ncid, 'g_basis', NF90_DOUBLE, [dimid_nloc, dimid_basisdim, dimid_enudim], varid) )
call check( nf90_put_var(out_ncid, varid, g_basis) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'g basis vectors') )
call check( nf90_put_att(out_ncid, varid, 'basis arrangement', '[east, north, up]') )

call check( nf90_def_var(out_ncid, 'ier', NF90_DOUBLE, [dimid_nloc], varid) )
call check( nf90_put_var(out_ncid, varid, ier_list) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'Apex error return. 0 good, 1 bad') )

! Close file
call check( nf90_close(out_ncid) )

!=============
contains
!=============
subroutine check(status)
integer, intent ( in) :: status

if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop "Stopped"
end if
end subroutine check  


end program matlab_apex_mall