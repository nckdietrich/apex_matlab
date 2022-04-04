function [qdlons, qdlats] = run_apex_mall(date, gdlons, gdlats, alts)
%{
Matlab function to run apex fortran code, specifically the apex_mall
subroutine which converts from apex cooredinates from geodetic coordinates.

Input:
    date  - Date in year and fraction (e.g., 2018.5) 
    gdlons - geodetic longitude [deg]
    gdlats - geodetic latitude [deg]
    alts  - altitude [km]

Output: 
    qdlons - Quasi-dipole longitude
    qdlats - Quasi-dipole latitude

Notes:
    Input a single date
    Ensure all longitudes, latitude and altitude arrays are equal length

Author: Nick Dietrich
Version: 4/1/22
%}

% Make column vector
gdlats = reshape(gdlats, [], 1);
gdlons = reshape(gdlons, [], 1);
alts = reshape(alts, [], 1);

nloc = length(gdlats);
% Check input vector lengths are consistent
if (length(gdlats) ~= length(gdlons)) && (length(gdlats) ~= length(alts))
    qdlons = [];
    qdlats = [];
    fprintf('Input vectors must be equal length \n');
    return
end

%% Write out to netcdf file
nfile_out = 'geodetic_pos.nc';
if exist(nfile_out,'file')
    delete(nfile_out);
end

% All positions
nccreate(nfile_out, 'loc',...
    'Dimensions', {'nloc',nloc, 'locDim', 3},...
    'Format','classic')
ncwrite(nfile_out, 'loc', [gdlons, gdlats, alts])
ncwriteatt(nfile_out,'loc','long_name','Locations in geodetic coordinates');
ncwriteatt(nfile_out,'loc', 'Arrangement', '[lon, lat, alt]')

% Date
nccreate(nfile_out, 'date',...
    'Dimensions', {'time',1},...
    'Format','classic')
ncwrite(nfile_out, 'date', date)
ncwriteatt(nfile_out,'date','long_name',' Date in year and fraction');

%% Running apex fortran function
func_path = what('apex_matlab');
system(strcat(func_path.path, '/matlab_apex_mall'));

%% Reading in fortran produced netcdf file
nfile_in = 'apex_pos.nc';
quasidipole_pos = ncread(nfile_in, 'qdipole_loc');
qdlons = quasidipole_pos(:,1);
qdlats = quasidipole_pos(:,2);

% Clean-up
delete(nfile_out);
% delete(nfile_in);
% cd(oriFolder);

end