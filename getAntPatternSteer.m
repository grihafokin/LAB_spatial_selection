% function for calculating the phase distribution of antenna array elements 
% input parameters:
%   antElPos -     array of coordinates for AA elements [x,y,z], m
%   f -            carrier frequency, Hz
%   azAng, elAng - direction (azimuth, elevation), in which the phase 
%            distribution of antenna array elements is calculated, degrees
function w = getAntPatternSteer(antElPos, f, azAng, elAng)
c = physconst('LightSpeed');
% vector, specifying the direction of arrival of the signal at the AA
incidentDir = [-cosd(elAng).*cosd(azAng);...
               -cosd(elAng).*sind(azAng);...
               -sind(elAng)];
% calculation of the projection of the direction vector onto the vectors, 
% connecting the origin of coordinates and position of the antenna array 
% elements, recalculation of this projection into a delay
tau = antElPos*incidentDir/c;
% calculation of phase distribution (in complex coefficient format)
% based on delay tau and carrier frequency f
w = exp(-1i*2*pi*f*tau);
end