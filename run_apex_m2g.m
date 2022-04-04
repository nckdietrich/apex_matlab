function [gdlons, gdlats] = run_apex_m2g(date, aplons, aplats, alts)
%{
Matlab function to run apex fortran code, specifically the apex_m2g
subroutine which converts modified apex coordinates to geodetic coordinates.

Input:
    date  - Date in year and fraction (e.g., 2018.5) 
    aplons - Quasi-dipole longitude
    aplats - Quasi-dipole latitude    
    alts  - altitude [km]

Output: 
    gdlons - geodetic longitude [deg]
    gdlats - geodetic latitude [deg]

Notes:
    Input a single date
    Ensure all longitudes, latitude and altitude arrays are equal length

Author: Nick Dietrich
Version: 4/4/22
%}

% Make column vector
aplats = reshape(aplats, [], 1);
aplons = reshape(aplons, [], 1);
alts = reshape(alts, [], 1);

nloc = length(aplats);
% Check input vector lengths are consistent
if (length(aplats) ~= length(aplons)) && (length(aplats) ~= length(alts))
    gdlons = [];
    gdlats = [];
    fprintf('Input vectors must be equal length \n');
    return
end

%% Write out to netcdf file
nfile_out = 'apex_pos.nc';
if exist(nfile_out,'file')
    delete(nfile_out);
end

% All positions
nccreate(nfile_out, 'loc',...
    'Dimensions', {'nloc',nloc, 'locDim', 3},...
    'Format','classic')
ncwrite(nfile_out, 'loc', [aplons, aplats, alts])
ncwriteatt(nfile_out,'loc','long_name','Locations in quasi-dipole coordinates');
ncwriteatt(nfile_out,'loc', 'Arrangement', '[lon, lat, alt]')

% Date
nccreate(nfile_out, 'date',...
    'Dimensions', {'time',1},...
    'Format','classic')
ncwrite(nfile_out, 'date', date)
ncwriteatt(nfile_out,'date','long_name',' Date in year and fraction');

%% Running apex fortran function
func_path = what('apex_matlab');
system(strcat(func_path.path, '/matlab_apex_m2g'));

%% Reading in fortran produced netcdf file
nfile_in = 'geodetic_pos.nc';
geodetic_pos = ncread(nfile_in, 'geodetic_loc');
gdlons = geodetic_pos(:,1);
gdlats = geodetic_pos(:,2);

% Clean-up
delete(nfile_out);
% delete(nfile_in);

end