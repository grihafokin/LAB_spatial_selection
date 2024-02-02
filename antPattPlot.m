% function for plotting antenna array (AA) radiation pattern (ARP)
% antElPos   - array of AA elements coordinates [x,y,z], m
% f          - carrier frequency, Hz
% nodeStruct - gNB or UE node structure with parameters
% angStep    - step of the grid of angles in azimuth and elevation
% scl        - scaling factor for visualization of antenna array pattern
% backLobe   - use backlobe suppression
% plt2D      - display ARP in 2D at an elevation angle slice equal to 0
% [x, y, z]  - surface coordinates for displaying ARP (dimensionless value,
%                                                   used for display only)
function [x, y, z, pp] = antPattPlot(antElPos, f, nodeStruct, ...
                                     angStep, scl, backLobe, plt2D)   
if (nargin == 6)
    plt2D = 0;
end
    
% grid of angles in azimuth and elevation for calculating the ARP
azA = 0:angStep:360;
if (plt2D == 0)
    elA = -90:angStep:90;
else
    elA = 0;
end
% % initialization of arrays for storing ARP values
azN = length(azA);
elN = length(elA);
x = zeros(elN, azN);
y = zeros(elN, azN);
z = zeros(elN, azN);
pp = zeros(elN, azN);
    
for i=1:elN % loop through an array of elevation angles
    for j=1:azN % loop through an array of azimuth angles     
        % calculation of ARP value for the i-th elevation and j-th azimuth
        % angle, account vector of steering coefficients (nodeStruct.Steer)
        if (plt2D == 0)
            p = getAntPatternG(antElPos, f, azA(j), elA(i),...
                nodeStruct.Steer, backLobe)*scl;
        else
            p = getAntPatternG(antElPos, f, azA(j),...
                nodeStruct.UEPolarCoord(2), ...
                nodeStruct.Steer, backLobe)*scl;
        end
        pp(i,j) = p;

        % recalculation of ARP values from polar to rectangular coordinates 
        % taking into account the orientation of AA (nodeStruct.AntOrient)
        r = p.*cosd(elA(i));
        if (plt2D == 0)
            xyz = nodeStruct.AntOrient*[r.*cosd(azA(j)); ...
                r.*sind(azA(j)); p.*sind(elA(i))];
        else
          xyz = [[cosd(-nodeStruct.AntDir(1)),sind(-nodeStruct.AntDir(1)); 
             -sind(-nodeStruct.AntDir(1)), cosd(-nodeStruct.AntDir(1))]*...
             [r.*cosd(azA(j)); r.*sind(azA(j))];0];
        end
        x(i,j) = xyz(1);
        y(i,j) = xyz(2);
        z(i,j) = xyz(3);
    end % for j=1:azN       % loop through an array of azimuth angles 
end % for i=1:elN           % loop through an array of elevation angles
end