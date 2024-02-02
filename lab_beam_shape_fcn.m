function [X,Y,Z] = lab_beam_shape_fcn(antType,Nel,stdCoords,...
                                        antPattCntrl,win_type,scl)
c = physconst('LightSpeed');
f = 30e9;       % carrier frequency in mmWave (millimetre wave) range, Hz
lamb = c/f;     % wavelength, m
da = 0.5*c/f;   % distance between antenna array (AA) elements, m
useAntUE = 0;   % use AA BF at the UE (antPattCntrl=0 mode only)
backLobe = 1;   % suppression of ARP back lobe (logical only for URA)
% period of the beamforming (BF)  procedure, s
Ta = 0; % when Ta = 0, beam is formed at each calculation step
antElPos = createAnt(antType, Nel, da); % construction of AA
NelFull = size(antElPos, 1);            % total number of AA elements
% scenario selection:
% 1- gNBs are located on the same side relative to the UE's motion path
% 2- gNBs are located on different sides relative to the UE's movement path
[gNB, ueNode, d, T, v] = createScenarion(1);
Nd = length(d);

Nnb = length(gNB);    % gNB number
Nue = length(ueNode); % UE number
% number of calculation points=number of coordinate points of UE trajectory
N = length(ueNode(1).Trajectory(:,1));
trajArray = [ueNode.Trajectory];       % UE coordinate array [N x 3*Nue]
gNBcoords = [gNB(:).Coords].';         % gNB coordinate array [Nnb x 3]
% initializing arrays to store angles of departure from each gNB
% to each UE (needed to eliminate the error of reading a non-existent
% array with some model settings)
azAng = zeros(Nue, Nnb);
elAng = zeros(Nue, Nnb);
azAngUE = zeros(Nue, Nnb);
elAngUE = zeros(Nue, Nnb);
% CYCLE BY THE NUMBER OF CALCULATION POINTS
for i=1:N % cycle by number of calculation points
    % array of coordinates of all UEs for the i-th calculation point
    ueCoordsi = reshape(trajArray(i, :).', 3, Nue).';
    % introducing an error into the estimation 
    % of UE coordinates according to stdCoords
    ueCoordsiErr = ueCoordsi;
    ueCoordsiErr(:,1:2) = ueCoordsiErr(:,1:2) + stdCoords*randn(2,2);
    if (mod(i, round(Ta/T)) == 1 || Ta == 0)
        % initializing arrays for storing 
        % angles of departure (AOD) from each gNB to each UE
        azAng = zeros(Nue, Nnb);  % azimuth angle of departure (AOD)
        elAng = zeros(Nue, Nnb);  % elevation angle of departure (AOD)
        azAngT = zeros(Nue, Nnb); % true azimuth angle of departure (AOD)
        elAngT = zeros(Nue, Nnb); % true elevation angle of departure (AOD)
        % initialization of arrays for storing angles of arrival (AOA)
        % for each UE from each gNB (used if there is an AA on the UE)
        azAngUE = zeros(Nue, Nnb);
        elAngUE = zeros(Nue, Nnb);
        for j=1:Nue % cycle by number of UEs
            % vector specifying the direction from gNB to UE 
            % in the global coordinate system x,y,z
            diffCoord = ueCoordsiErr(j,:) - gNBcoords;
            diffCoordT = ueCoordsi(j,:) - gNBcoords;
            for n=1:Nnb % cycle by gNB number
                % vector specifying the direction from gNB to UE in local
                % coordinate system of the gNB antenna array (i.e., taking 
                % into account the position of the gNB antenna array)
                dirVect = gNB(n).AntOrient.'*diffCoord(n,:).';
                dirVectT = gNB(n).AntOrient.'*diffCoordT(n,:).';
                % calculation of the angles of departure (AOD) 
                % from the n-th gNB to the j-th UE
                azAng(j, n) = rad2deg(atan2(dirVect(2), dirVect(1)));
                elAng(j, n) = rad2deg(atan2(dirVect(3), ...
                    sqrt(sum(dirVect(1:2).^2))));
                azAngT(j, n) = rad2deg(atan2(dirVectT(2), dirVectT(1)));
                elAngT(j, n) = rad2deg(atan2(dirVectT(3), ...
                    sqrt(sum(dirVectT(1:2).^2))));
                % calculation of the angles of departure (AOD)
                % for UE antenna array
                if ( useAntUE == 1)
                    % vector specifying the direction from the UE to gNB 
                    % in local coordinate system of UE antenna array, i.e. 
                    % taking into account position of the UE antenna array
                    dirVect = -ueNode(j).AntOrient.'*diffCoord(n,:).';
                    % calculation of the angles of departure (AOD) 
                    % from the j-th UE to the n-th gNB
                    azAngUE(j, n) = rad2deg(atan2(dirVect(2), dirVect(1)));
                    elAngUE(j, n) = rad2deg(atan2(dirVect(3), ...
                        sqrt(sum(dirVect(1:2).^2))));
                end
            end
        end
        % for an adaptive scheme, because it doesn't use information
        % about UE location
        if (antPattCntrl == 2)
            azAng = azAngT;
            elAng = elAngT;
        end
        for n=1:Nnb % cycle by gNB number
            cVect = zeros(NelFull, Nue);
            gVect = zeros(1, Nue);
            for j=1:Nue % cycle by UE number      
                % calculation of the steering vector of phase distribution 
                % of the n-th gNB antenna array (AA) to the j-th UE
                cVect(:,j) = ...
                 getAntPatternSteer(antElPos, f, azAng(j, n), elAng(j, n));
                % storage of data about the served UE (serving gNB number
                % is specified in the servgNB parameter for each UE)
                if (n == ueNode(j).servgNB)
                    % calculation of the vector of steering coefficients
                    % of the n-th gNB antenna array, serving the j-th UE
                    gNB(n).Steer = getAntPatternSteer(antElPos, f, ...
                        azAng(j, n), elAng(j, n))/NelFull;
                    % entry 1 in vector g by index of the j-th UE 
                    % to set the maximum in this direction
                    gVect(j) = 1;
                    % storing UE location information in polar
                    % coordinates (azimuth, elevation, distance)
                    gNB(n).UEPolarCoord = [azAng(j, n), elAng(j, n),...
                    sqrt(sum((ueCoordsiErr(j,:).' - gNB(n).Coords).^2))];
                end
                % calculation of the vector of steering coefficients
                % of the j-th UE antenna array, working with the n-th gNB
                if (n == ueNode(j).servgNB && useAntUE == 1)
                    ueNode(j).Steer = getAntPatternSteer(antElPos, f, ...
                        azAngUE(j,n), elAngUE(j,n));
                end
            end
            % calculation of the steering vector of coefficients of n-th
            % gNB antenna array, serving j-th UE with a beam control scheme 
            % different from antPattCntrl=0
            switch antPattCntrl
                % LCMV
                case 1
                    % calculation of the vector of coefficients using the 
                    % LCMV algorithm for control maximum and null beam ARP
                    gNB(n).Steer = cVect*pinv(cVect'*cVect)*gVect';
                    % normalization of coefficients
                    gNB(n).Steer = ...
                        gNB(n).Steer/max(abs(gNB(n).Steer))/NelFull;
                % SMI
                case 2
                    K = 100; % preamble length
                    % steering vector of antenna array for the SOI signal
                    sv_s = cVect(:,[ueNode.servgNB] == n);
                    % steering vector of antenna array for the SNOI signal
                    sv_j = cVect(:,~([ueNode.servgNB] == n));
                    % simulation of signal and interference reception 
                    % from specified SOI and SNOI directions
                    rx_sum = ones(K,1)*sv_s.' + ...
                        (randn(K,1) + 1i*randn(K,1))/sqrt(2)*sv_j.';
                    % correlation matrix of signals from AA elements
                    Rxx = rx_sum'*rx_sum;
                    % correlation vector of AA element signals
                    % with reference preamble
                    rxd = rx_sum'*ones(K,1);
                    % calculation of vector of coefficients according to 
                    % the SMI algorithm to suppress interference
                    gNB(n).Steer = conj(pinv(Rxx)*rxd);
                % beamshaping
                case 3
                    % ARP beam width multiplier
                    switch win_type
                        case 0
                            scl = 0.45;
                        case 1
                            scl = 0.1;
                        case 2
                            scl = 2.2;
                    end                 
                    % ARP beam width, calculated from the standard 
                    % deviation of the UE coordinate estimate
                    BW = 2*atan2d(stdCoords,gNB(n).UEPolarCoord(3));
                    % calculation of the AA vector of weight coefficients 
                    w = beamshapingWeight(win_type, BW, ...
                        gNB(n).UEPolarCoord(1), Nel, scl);
                    % apply coefficient w, if the required beam width is:
                    % not less than the minimum theoretical beam width,
                    % equal to 0.891*lamb/Nel/da; 
                    % amplitude of the sum of coefficients>0.1
                    if ~(BW*scl <= rad2deg(0.891*lamb/Nel/da) ||...
                            any(isnan(w)) || sum(w) < 0.1)
                        if (antType == 1)
                            % calculation of the vector of weighting 
                            % coefficients of URA for vertical AA elements
                            antElPosElev = [zeros(Nel,1), ...
                                antElPos(1:Nel,2), zeros(Nel,1)];
                            wEl = getAntPatternSteer(antElPosElev, ...
                                f, gNB(n).UEPolarCoord(2), 0)/Nel;
                            % vector of weighting coefficients of the URA,
                            % obtained by multiplying horizontal and
                            % vertical AA elements coefficients
                            w = w*wEl.'/Nel;
                            w = w(:);
                        end                        
                        gNB(n).Steer = w;
                    end
            end % switch antPattCntrl
        end % for n=1:Nnb     
    end % if (mod(i, round(Ta/T)) == 1 || Ta == 0)
    % array for temporary storage of values of
    % received power per UE from each gNB
    gNBpwr = zeros(Nnb, 1);
    % calculation of received power ratio from the serving gNB
    % to power, received from neighboring (interfering) gNB for each UE
    for j=1:Nue % cycle by UE number
        for n=1:Nnb % cycle by gNB number
            % calculation of the gain of the UE receiving antenna array
            if ( useAntUE == 1 && j == 1)
                % gain in the presence of antenna array on the UE, 
                % the beam of which is directed to the serving gNB
                gUE = getAntPatternG(antElPos, f, ...
                    azAngUE(n), elAngUE(n), ueNode(1).Steer, backLobe).^2;
            else
                % gain without antenna array on the UE
                gUE = 1;
            end
            % calculation of the power, received at j-th UE from n-th gNB, 
            % taking into account beamforming on UE & gNB,excluding range
            gNBpwr(n) = gUE*getAntPatternG(antElPos, f, ...
                azAngT(j, n), elAngT(j, n), gNB(n).Steer, backLobe).^2;
        end
        % calculating the distance from the j-th UE to each gNB
        diffCoord = ueCoordsi(j,:) - gNBcoords;
        distSpace = sqrt(sum(diffCoord.^2,2));
        % calculation of the power, received at the j-th UE from each gNB, 
        % taking into account range (losses are calculated using the free 
        % space path loss attenuation model)
        gNBpwr = pow2db(gNBpwr) - fspl(distSpace,c/f);
        % calculation of received power ratio from serving gNB
        % to the power, received from the neighboring gNB for the j-th UE
        ueNode(j).SNR(i) = gNBpwr(ueNode(j).servgNB) - ...
            sum(gNBpwr(1:end ~= ueNode(j).servgNB));
    end
end % for i=1:N 
% preparing arrays of coordinates and array of values 
% of signal-to-interference ratio of 1st UE for display
X = reshape(ueNode(1).Trajectory(:,1), [], Nd).';
Y = reshape(ueNode(1).Trajectory(:,2), [], Nd).';
Z = reshape(ueNode(1).SNR, [], Nd).';
end