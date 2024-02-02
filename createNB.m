% gNB structure creation function
% input parameters:
%   coords - gNB coordinates [x,y,z], m
%   antDir - gNB antenna array orientation (direction of symmetry axis) 
%            [azimuth, tilt], degrees
function gNB = createNB(coords, antDir)
gNB.Coords = coords.'; % gNB coordinates
gNB.AntDir = antDir.'; % gNB antenna array (FF) orientation
% rotation matrix according to antDir values, used for recalculation of 
% direction vectors from gNB global coordinates to gNB AA local coordinates
gNB.AntOrient = rotz(antDir(1))*roty(-antDir(2));
gNB.Steer = 1;
gNB.UEPolarCoord = [0, 0, 0];
end