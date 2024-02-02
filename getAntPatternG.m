% function for calculation of antenna array (AA) gain 
% input parameters:
%   antElPos -     array of AA elements coordinates [x,y,z], m
%   f -            carrier frequency, Hz
%   azAng, elAng - direction (azimuth, elevation),  
%                  in which the gain is calculated, degrees
%   wS -           vector of steering coefficients
%   backLobe -     use backlobe suppression
% output parameters:
%   g - AA gain in absolute values in the direction [azAng, elAng]
function g = getAntPatternG(antElPos, f, azAng, elAng, wS, backLobe)
% beamforming, i.e. application of wS to coefficients, describing
% phase distribution of AA elements for a given direction [azAng, elAng]
% w = wS'*getAntPatternSteer(antElPos, f, azAng, elAng);
w = sum(conj(wS).*getAntPatternSteer(antElPos, f, azAng, elAng));
g = abs(w); % AA gain in amplitude terms
% backlobe suppression via reduction of AA gain
% in the angle sector from 90 to 270 degrees
if (backLobe == 1)
    if (azAng > 90 && azAng < 270)
        azAng = wrapTo360(azAng) - 90;
        % as coefficient, suppressing AA pattern in the angle sector from 
        % 90 to 270 degrees we use ellipse radius with eccentricity equal 
        % to 0.9999; the closer the angle to 180, the higher is suppression
        g = g*sqrt(cosd(azAng)^2 + (sqrt(1-0.9999^2)*sind(azAng))^2);
    end
end
end