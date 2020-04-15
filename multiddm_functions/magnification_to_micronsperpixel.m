function [microns_per_pixel, mag_str] = magnification_to_micronsperpixel(video_fname, json_filename)
%magnification_to_micronsperpixel returns um/px given string with magnification
%   Parse the video_fname looking for a string that is recognisable as
%   magnification string (e.g. 20X_), read the corresponding value from
%   json_filename, return the correct value for um/px (or px2mum in
%   DDM_Analysis convention)
%   This is not a very smart function, so I didn't write a very nice regexp
%   to deal with all possible cases. 
%   It scans for a number (with decimals) followed by 'X' (or 'x')
%   It can deal with a string like '60X_1.5X' like used some times when the 
%   additional 1.5x lens in a nikon microscope is engaged (magnification
%   90X is used in this case)

% text match
regex = '\d+\.?\d*?(?=X)';
match = regexpi(video_fname, regex, 'match', 'all');

% deal with how many matches
if numel(match) == 1
    mag = str2double(match{1});
elseif numel(match) == 2
    mag = str2double(match{1}) * str2double(match{2});
else
    warning(['Compatible magnification string not found in ', video_fname])
    microns_per_pixel = -1;
end % if

% convert back into a magnification string compatible with matlab
mag_str = [num2str(mag), 'X'];
mag_str_tofind = matlab.lang.makeValidName(mag_str);

% read json
known_mags = jsondecode(fileread(json_filename));

% check we have a match
is_known_mag = ismember(mag_str_tofind, fieldnames(known_mags));
if ~is_known_mag
    warning(['Couldn''t find a magnification match in ', json_filename])
    microns_per_pixel = -1;
else 
    % and read the found match
    microns_per_pixel = known_mags.(mag_str_tofind);
end % if

