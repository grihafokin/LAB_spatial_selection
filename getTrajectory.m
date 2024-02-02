% function for constructing the trajectory of the UE movement
% input parameters:
%   trajPnts - initial and final values of UE trajectory along X, Y, Z, m
%   v - UE movement speed, m/s
%   T - UE coordinate measurement period, s
function xyUN = getTrajectory(trajPnts, v, T)
% total time of the UE movement scenario
times = [0; sqrt(sum((trajPnts(2:end, :) - trajPnts(1:end-1, :)).^2,2))/v];
% time at each reference point (trajPnts) describing trajectory of movement
elapsedTime = zeros(1, length(times));    
for i=1:length(times)
    elapsedTime(i) = sum(times(1:i));
end
% built-in functions for generating coordinates
% of UE trajectory at specified reference points
ts = trackingScenario('UpdateRate', 1/T);
target = platform(ts);
traj = waypointTrajectory('Waypoints', trajPnts, ...
    'TimeOfArrival', elapsedTime, 'SampleRate', 1/T);
target.Trajectory = traj;
r = record(ts);
posUN = [r(:).Poses];
% coordinates of UE trajectory points
xyUN = vertcat(posUN.Position);
end