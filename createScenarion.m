% function for creating a scenario: gNB position, UE movement trajectories
function [gNB, ueNode, d, T, v] = createScenarion(sceneN, varargin)
T = 0.1;  % measurement period
v = 10;   % UE device movement speed, m/s
switch sceneN
    % gNBs are located on the same side relative to the trajectory
    % of UE movement; UEs move in parallel at a distance d
    case 1 
        % distance between the motion trajectories of two UEs;
        % in this scenario the UE trajectories are spaced by d in Y
        if (nargin == 2)
            d = varargin{1};
        else
            d = 0;
        end
        % creating array of two gNB structures; gNB parameters in createNB
        gNB(1) = createNB([25, 50, 5], [-90, -1]);
        gNB(2) = createNB([125, 50, 5], [-90, -1]);
        %creating array of two UE structures; UE parameters in createUEnode
        ueNode(1) = createUEnode([0; 150], [0; 0], 0, 1, v, T, [90, 0]);
        ueNode(2) = createUEnode([150; 0], [0; 0], 0, 2, v, T, [90, 0]); 
        % constructing a set of UE motion trajectories for different values 
        % of d;trajectory UE2 does not change, trajectory UE1 shifts by d
        Nd = length(d); % number of d values
        trajUE1 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        trajUE2 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        for i=1:Nd
            trajUE1(:,:,i) = ueNode(1).Trajectory + [0, d(i), 0];
            trajUE2(:,:,i) = ueNode(2).Trajectory;
        end
        % arrangement of a set of UE trajectories into one common array of 
        % coordinates for each UE
        trajUE1 = permute(trajUE1, [1 3 2]);
        trajUE1 = reshape(trajUE1, [],3 ,1);
        trajUE2 = permute(trajUE2, [1 3 2]);
        trajUE2 = reshape(trajUE2, [],3 ,1);
        % saving trajectory in UE structure
        ueNode(1).Trajectory = trajUE1;
        ueNode(2).Trajectory = trajUE2;

    % gNBs are located on different sides relative to the trajectory
    % of UE movement; UEs move in parallel at a distance d
    case 2      
        % distance between the motion trajectories of two UEs;
        % in this scenario the UE trajectories are spaced by d in Y
        d = 0:1:10;
        % creating an array of two gNB structures
        gNB(1) = createNB([25, 50, 5], [-90, -1]);
        gNB(2) = createNB([125, 0, 5], [90, -1]);
        % creating an array of two UE structures
        ueNode(1) = createUEnode([0; 150], [20; 20], 0, 1, v, T, [90, 0]);
        ueNode(2) = createUEnode([150; 0], [20; 20], 0, 2, v, T, [90, 0]);
        % constructing a set of UE motion trajectories for different values 
        % of d; trajectory UE2 does not change, trajectory UE1 shifts by d
        Nd = length(d); % number of d values
        trajUE1 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        trajUE2 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        for i=1:Nd
            trajUE1(:,:,i) = ueNode(1).Trajectory + [0, d(i), 0];
            trajUE2(:,:,i) = ueNode(2).Trajectory;
        end
        % arrangement of a set of UE trajectories into one common array of 
        % coordinates for each UE
        trajUE1 = permute(trajUE1, [1 3 2]);
        trajUE1 = reshape(trajUE1, [],3 ,1);
        trajUE2 = permute(trajUE2, [1 3 2]);
        trajUE2 = reshape(trajUE2, [],3 ,1);
        % saving trajectory in UE structure
        ueNode(1).Trajectory = trajUE1;
        ueNode(2).Trajectory = trajUE2;
   
    % gNBs are located on the same side relative to the trajectory
    % of UE movement; UEs move one after another at a distance d
    case 3
        % distance between the motion trajectories of two UEs;
        % in this scenario the UE trajectories are spaced by d in X
        d = 0:2:150;
        % creating an array of two gNB structures
        gNB(1) = createNB([25, 50, 5], [-90, -1]);
        gNB(2) = createNB([125, 50, 5], [-90, -1]);
        % creating an array of two UE structures
        ueNode(1) = createUEnode([0; 150], [0; 0], 0, 1, v, T, [90, 0]);
        ueNode(2) = createUEnode([0; 150], [0; 0], 0, 2, v, T, [90, 0]);       
        % constructing a set of UE motion trajectories for different values
        % of d; trajectory UE2 does not change, trajectory UE1 shifts by d
        Nd = length(d); % число значений d
        trajUE1 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        trajUE2 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        for i=1:Nd
            trajUE1(:,:,i) = ueNode(1).Trajectory + [d(i), (i-1)/1e10, 0];
            trajUE2(:,:,i) = ueNode(2).Trajectory;
        end
        % arrangement of a set of trajectories into one common array of 
        % coordinates for each UE
        trajUE1 = permute(trajUE1, [1 3 2]);
        trajUE1 = reshape(trajUE1, [],3 ,1);
        trajUE2 = permute(trajUE2, [1 3 2]);
        trajUE2 = reshape(trajUE2, [],3 ,1);
        % saving trajectory in UE structure
        ueNode(1).Trajectory = trajUE1;
        ueNode(2).Trajectory = trajUE2;
    otherwise
end
end