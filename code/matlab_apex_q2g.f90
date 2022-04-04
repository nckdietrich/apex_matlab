! Running apex_mall - Converts quasi-dipole coords to geodetic coords
! Author: Nick Dietrich
! Version 4/1/22

program matlab_apex_q2g

use apex
use netcdf

implicit none

!! Declaring 
character(len=128) :: file_in = 'quasidipole_pos.nc'
character(len=128) :: file_out = 'geodetic_pos.nc'
integer :: in_ncid, dimid, varid
integer :: out_ncid, dimid_nloc, dimid_locdim, dimid_basisdim, dimid_enudim
integer :: nloc, loc_ind, ier
real :: date, altmax, hr
real, allocatable :: qdipole_loc(:,:), qdlons(:), qdlats(:), alts(:)
! Apex output variables
real :: gdlat, gdlon
real, allocatable :: gd_loc(:,:), ier_list(:)

!! Running
print *, 'Converting Quasi-Dipole Coords to Geodetic Coords'
! Read in information
! Get vector length
call check( nf90_open(file_in, NF90_NOWRITE, in_ncid) )
call check( nf90_inq_dimid(in_ncid, 'nloc', dimid))
call check( nf90_inquire_dimension(in_ncid, dimid, len=nloc))
allocate(qdipole_loc(nloc,3))
! Locations [lon, lat, alt]
call check( nf90_inq_varid(in_ncid, 'loc', varid) )
call check( nf90_get_var(in_ncid, varid, values=qdipole_loc) )
! date
call check( nf90_inq_varid(in_ncid, 'date', varid) )
call check( nf90_get_var(in_ncid, varid ,values=date))
! close
call check( nf90_close(in_ncid) )

! Running apex code
hr = 90 ! Typical value used in TIEGCM [km]
allocate(qdlons(nloc))
allocate(qdlats(nloc))
allocate(alts(nloc))
qdlons = qdipole_loc(:,1)
qdlats = qdipole_loc(:,2)
alts = qdipole_loc(:,3)
altmax = MAXVAL(alts)

! Set up interp arrays (takes longest to run) + allocate
call apex_setup(date, altmax)
allocate(gd_loc(nloc,3))
allocate(ier_list(nloc))
! Run apex
print *, 'Looping through locations...'
do loc_ind = 1, nloc
    call apex_q2g( qdlats(loc_ind), qdlons(loc_ind), alts(loc_ind), &
        gdlat, gdlon, ier )

    ! Save to variables
    gd_loc(loc_ind,:) = [gdlon, gdlat, alts(loc_ind)]
    ier_list(loc_ind) = ier    
enddo

! Write out data (saving everyting)
print *, 'Writing out to file...'
call check( nf90_create(file_out, NF90_NETCDF4, out_ncid) )
! Define dimensions
call check( nf90_def_dim(out_ncid, 'nloc', nloc, dimid_nloc) )
call check( nf90_def_dim(out_ncid, 'locDim', 3, dimid_locdim))

! Writing out variables
call check( nf90_def_var(out_ncid, 'geodetic_loc', NF90_DOUBLE, [dimid_nloc, dimid_locdim], varid) )
call check( nf90_put_var(out_ncid, varid, gd_loc) )
call check( nf90_put_att(out_ncid, varid, 'long_name', 'Geodetic Corrdinates') )
call check( nf90_put_att(out_ncid, varid, 'arrangement', '[lon, lat, alt]') )
call check( nf90_put_att(out_ncid, varid, 'units', '[deg, deg, km]') )

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


end program matlab_apex_q2g