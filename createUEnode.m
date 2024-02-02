% UE structure creation function
% input parameters:
%   xPnts - initial and final values of the UE trajectory along X, m
%   yPnts - initial and final values of the UE trajectory along Y, m
%   zPnts - height UE (constant throughout the entire trajectory), m
%   servgNB - number of serving gNB (index in the array of gNB structures)
%   v -       UE movement speed, m/s
%   T -       UE coordinate measurement period, s
%   antDir -  UE antenna array orientation (symmetry axis direction)
%             [azimuth, inclination], degrees
function ue = createUEnode(xPnts, yPnts, zPnts, servgNB, v, T, antDir)
ue.servgNB = servgNB; % serving gNB number  
% creating an array of UE motion trajectory coordinates
trajPnts = [xPnts, yPnts];
trajPnts(:,3) = zPnts;
ue.Trajectory = getTrajectory(trajPnts, v, T);    
% initializing an array of signal-to-interference ratio values,
% is calculated at each point of the UE trajectory
ue.SNR = zeros(length(ue.Trajectory(:,1)), 1);    
% rotation matrix according to antDir values, used for recalculation 
% direction vectors from global coordinates to the UE antenna array 
% coordinate system
ue.AntOrient = rotz(antDir(1))*roty(-antDir(2));
% вектор направляющих коэфф. АР (поумолчанию АР не направлена)
% vector of antenna array steering coefficients; 
% by default AA is not directed
ue.Steer = 1;
end