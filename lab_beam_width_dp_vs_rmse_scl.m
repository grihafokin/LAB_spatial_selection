% 
clear all; close all; clc;
stdCoordsArr = [1, 3, 5, 10];  % RMSE of UE coordinate estimates, m
NellArr = 4:1:20;              % number of AA elements in one dimension
% selection of antenna array type
% 1 - planar or uniform rectangural antenna array (URA)
% 2 - uniform linear antenna array (ULA)
% 3 - uniform circular antenna array (UCA); not supported
antType = 1;
% select ARP beam shape for beam width control algorithm (antPattCntrl = 3)
% 0 - Gauss window
% 1 - raised cosine window
% 2 - rectangular window
win_typeArr = [2, 0];
antTypeCmt = ["URA", "ULA", "UCA"];
for ww=1:length(win_typeArr)
    win_type = win_typeArr(ww);
    switch win_type
        case 0
            sclArr = 0.1:0.1:3.0; % Gauss window
        case 1
            sclArr = 0.0:0.1:0.9; % Raised cosine window;
        case 2
            sclArr = 0.1:0.1:3.0; % Rectangular window
    end
    c = physconst('LightSpeed');
    f = 30e9;       % carrier frequency, Hz
    lamb = c/f;     % wavelength, m
    da = 0.5*c/f;   % distance between antenna array (AA) elements, m
    Z = [];
    for s=stdCoordsArr
        Zn = [];
        for nel=NellArr
            Zsc = [];
            for sc=sclArr
                rng('default');
                stdCoords = s; % RMSE of UE coordinate estimates, m
                Nel = nel;     % number of AA elements in one dimension
                antElPos=createAnt(antType, Nel, da); % construction of AA
                NelFull = size(antElPos, 1);  % total number of AA elements
                N = 100;  % number of calculation points
                HgNB = 0; % gNB antenna array height
                Due = 50; % gNB to UE 2D distance
                % gNB parameter structure
                gNB = createNB([0, 0, HgNB], [0, 0]);
                gNB.Steer = zeros(NelFull, 2);
                ueRxPwr = zeros(2, N);
                ueCoord = [Due, 0, 0];         % UE coordinates
                gNBcoords = [gNB(:).Coords].'; % array of gNB coordinates
                % gNB-UE 3D distance
                distSpaceT = sqrt(sum((gNB.Coords-ueCoord.').^2)); 
                BW = 2*atan2d(stdCoords,distSpaceT);% ARP HPBW
                stAng = 0;                          % ARP maximum direction
                % vector specifying the direction from gNB to UE 
                % in the global coordinate system x,y,z
                diffCoord = ueCoord - gNBcoords;
                % vector specifying the direction from gNB to UE
                % in the local coordinate system of the gNB antenna array, 
                % i.e., taking into account the position of the gNB AA
                dirVect = gNB.AntOrient.'*diffCoord.';
                % calculate angle of departure (AOD) from gNB to UE 
                azAng = rad2deg(atan2(dirVect(2), dirVect(1)));
                elAng = rad2deg(atan2(dirVect(3), ...
                    sqrt(sum(dirVect(1:2).^2))));
                % calculation of vector of steering coefficients of gNB AA
                gNB.Steer(:,1) = getAntPatternSteer(antElPos, f,...
                    azAng, elAng)/NelFull;
                % calculation of the vector of weight coefficients of AA
                [w, azAngP, antPattP] = beamshapingWeight(win_type, BW,...
                                            stAng, Nel, sc);
                % apply w coefficients if the required beam width is not 
                % less, than minimum theoretical width 0.891*lamb/Nel/da
                if (BW*sc <= rad2deg(0.891*lamb/Nel/da)...
                        || any(isnan(w)) || sum(w) == 0)
                    w = gNB.Steer(:,1);
                else
                    if (antType == 1)
                        % calculation of the vector of weight coefficients 
                        % for vertical elements of URA
                        w = repmat(w, Nel, 1)/Nel;
                        w = w(:);
                    end
                end
                gNB.Steer(:,2) = w;
                for i=1:N % cycle by number of calculation points
                    % introducing an error into the estimation 
                    % of UE coordinates according to stdCoords
                    ueCoordErr = ueCoord;
                    ueCoordErr(1:2) = ueCoordErr(1:2) +...
                        stdCoords*randn(size(ueCoord(1:2)));
                    % vector specifying the direction from gNB to UE in 
                    % the global coordinate system, accounting RMSE of 
                    % UE coordinate estimate
                    diffCoordT = ueCoordErr - gNBcoords;
                    % vector specifying direction from gNB to UE in local
                    % coordinate system of gNB antenna array (i.e., taking 
                    % into account the position of the gNB antenna array)
                    dirVectT = gNB.AntOrient.'*diffCoordT.';
                    % calculation of the AOD from gNB to UE
                    azAngT = rad2deg(atan2(dirVectT(2), dirVectT(1)));
                    elAngT = rad2deg(atan2(dirVectT(3), ...
                        sqrt(sum(dirVectT(1:2).^2))));
                    % calculation of received power from the serving gNB, 
                    % accounting BF at gNB and w/o distance accounting 
                    gNBpwr = [getAntPatternG(antElPos, f, ...
                        azAngT, elAngT, gNB.Steer(:,1), 0).^2;...
                        getAntPatternG(antElPos, f, ...
                        azAngT, elAngT, gNB.Steer(:,2), 0).^2];
                    % calculation of distance from UE to gNB
                    diffCoord = ueCoordErr - gNBcoords;
                    distSpace = sqrt(sum(diffCoord.^2,2));
                    % calculation of the power, received at UE from gNB, 
                    % taking into account range (losses are calculated 
                    % using the free space path loss attenuation model)
                    gNBpwr(isnan(gNBpwr)) = gNBpwr(1);
                    gNBpwr = pow2db(gNBpwr) - fspl(distSpace,c/f);
                    ueRxPwr(:, i) = gNBpwr;
                end % for i=1:N 
                Zi = mean(ueRxPwr(2,:) - ueRxPwr(1,:));
                fprintf('rmse = %5.1f, nel=%3i, scl=%4.2f\n', s, nel, sc);
                Zsc = [Zsc, Zi];
            end
            Zn = [Zn;Zsc];
        end
        Z{end+1} = Zn;
    end
    figure;
    [X,Y] = meshgrid(NellArr, sclArr);
    for kk=1:size(Z,2)
        subplot(2,2,kk)
        if stdCoordsArr(kk)==1
            [~,hh]=contourf(X, Y, smoothdata(Z{kk}.')); 
        else
            [~,hh]=contourf(X, Y, smoothdata(Z{kk}.'),'ShowText','on'); 
        end
        c1 = colorbar; c1.Label.String = '\DeltaP, dB'; hold on;
        hh.LabelSpacing = 250;
        grid on; xlabel('N'); ylabel('s'); view([0, 90]); axis tight;
        title(sprintf('RMSE = %d m',stdCoordsArr(kk)));
    end
    if win_type == 0 
        sgtitle(sprintf('Gaussian ARP; %s', antTypeCmt(antType)));
    elseif win_type == 2
        sgtitle(sprintf('Rectangular ARP; %s', antTypeCmt(antType)));
    end   
end