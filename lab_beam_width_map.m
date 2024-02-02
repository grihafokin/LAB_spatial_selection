clear all; close all; clc;
% selection of antenna array type
% 1 - planar or uniform rectangural antenna array (URA)
% 2 - uniform linear antenna array (ULA)
% 3 - uniform circular antenna array (UCA); is not supported
antType = 1;
Nel = 20;          % number of AA elements in one dimension
c = physconst('LightSpeed');
f = 30e9;       % carrier frequency, Hz
lamb = c/f;     % wavelength, m
da = 0.5*c/f;   % distance between antenna array (AA) elements, m
antElPos = createAnt(antType, Nel, da); % construction of AA
NelFull = size(antElPos, 1);            % total number of AA elements
% select ARP beam shape for beam width control algorithm (antPattCntrl = 3)
% 0 - Gauss window
% 1 - raised cosine window
% 2 - rectangular window
win_typeArr = [2, 0];
antTypeCmt = ["URA", "ULA", "UCA"];
ueRxPwrPlt = [];
figNumber = 1;
stdCoordsArr=[10];
for s=stdCoordsArr
    subFigNumber = 1;
    for ww=win_typeArr
        stdCoords = s; % RMSE of UE coordinate estimates, m
        c = physconst('LightSpeed');
        f = 30e9;       % carrier frequency, Hz
        lamb = c/f;     % wavelength, m
        da = 0.5*c/f;   % distance between antenna array (AA) elements, m
        stdCoordsArr = 10; % RMSE of UE coordinate estimates along, m
        N = 100;  % number of calculation points
        HgNB = 0; % gNB antenna array height
        Due = 50; % gNB to UE distance
        % gNB parameter structure
        gNB = createNB([0, 0, HgNB], [0, 0]);
        gNB.Steer = zeros(NelFull, 2);
        ueRxPwr = zeros(2, N);
        ueCoord = [Due, 0, 0];         % UE coordinates
        gNBcoords = [gNB(:).Coords].'; % array of gNB coordinates
        distSpaceT = sqrt(sum((gNB.Coords-ueCoord.').^2)); % gNB-UE 3D dist
        BW = 2*atan2d(stdCoords,distSpaceT);        % ARP beam width
        stAng = 0;                                  % ARP maximum direction
        % vector specifying the direction from gNB to UE 
        % in the global coordinate system x,y,z
        diffCoord = ueCoord - gNBcoords;
        % vector specifying the direction from gNB to UE
        % in the local coordinate system of the gNB antenna array, 
        % i.e., taking into account the position of the gNB antenna array
        dirVect = gNB.AntOrient.'*diffCoord.';
        % calculate angle of departure (AOD) from gNB to UE 
        azAng = rad2deg(atan2(dirVect(2), dirVect(1)));
        elAng = rad2deg(atan2(dirVect(3), sqrt(sum(dirVect(1:2).^2))));
        % calculation of the vector of steering coefficients of gNB AA
        gNB.Steer(:,1)=getAntPatternSteer(antElPos,f,azAng,elAng)/NelFull;
        scl = 1.5; % scaling multiplier of ARP beam width
        % calculation of the vector of weight coefficients of antenna array
        [w, azAngP, antPattP] = beamshapingWeight(ww, BW, stAng, Nel, scl);
        % apply w coefficients if the required beam width is not less,
        % than the minimum theoretical width 0.891*lamb/Nel/da
        if (BW*scl <= rad2deg(0.891*lamb/Nel/da) || any(isnan(w)))
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
        rng('default');
        for i=1:N % cycle by number of calculation points
            % introducing an error into the estimation 
            % of UE coordinates according to stdCoords
            ueCoordErr = ueCoord;
            ueCoordErr(1:2) = ueCoordErr(1:2) + ...
                stdCoords*randn(size(ueCoord(1:2)));
            % vector specifying the direction from gNB to UE in the global 
            % coordinate system, accounting RMSE of UE coordinate estimate
            diffCoordT = ueCoordErr - gNBcoords;
            % vector specifying the direction from gNB to UE in local
            % coordinate system of the gNB antenna array (i.e., taking 
            % into account the position of the gNB antenna array)
            dirVectT = gNB.AntOrient.'*diffCoordT.';
            % calculation of the angles of departure (AOD) from gNB to UE
            azAngT = rad2deg(atan2(dirVectT(2), dirVectT(1)));
            elAngT=rad2deg(atan2(dirVectT(3),sqrt(sum(dirVectT(1:2).^2))));
            % calculation of received power from the serving gNB, 
            % accounting beamforming at gNB and w/o distance accounting 
            gNBpwr = [getAntPatternG(antElPos, f, ...
                        azAngT, elAngT, gNB.Steer(:,1), 0).^2;...
                        getAntPatternG(antElPos, f, ...
                        azAngT, elAngT, gNB.Steer(:,2), 0).^2];
            % calculation of distance from UE to gNB
            diffCoord = ueCoordErr - gNBcoords;
            distSpace = sqrt(sum(diffCoord.^2,2));
            % calculation of the power, received at the UE from gNB, 
            % taking into account range (losses are calculated using 
            % the free space path loss attenuation model)
            gNBpwr(isnan(gNBpwr)) = gNBpwr(1);
            gNBpwr = pow2db(gNBpwr) - fspl(distSpace,c/f);
            ueRxPwr(:, i) = gNBpwr;
        end % for i=1:N 
        % calculation of the ARP for the case without/with control of HPBW
        azAngPatt = -90:0.1:89;
        g = zeros(1, length(azAngPatt));
        gDef = zeros(1, length(azAngPatt));
        for i=1:length(azAngPatt)
            g(i) = getAntPatternG(antElPos, f, azAngPatt(i), 0, w, 0);
            gDef(i) = getAntPatternG(antElPos, f, ...
            azAngPatt(i), 0, ones(size(w))/NelFull, 0);
        end
        gNorm = g/max(g);
        % calculation of ARP half power beam width (HPBW) at -3 dB level
        [~,ind] = min(abs(gNorm - 1/sqrt(2)));
        hpbw = 2*abs(azAngPatt(ind));
        % calculation of coordinates of ARP curves for display,
        % rotated 90 degrees for better visualization
        alphPatt = azAngPatt + 90;
        xPatt = cosd(alphPatt).*gNorm*(distSpaceT);
        yPatt = sind(alphPatt).*gNorm*(distSpaceT);
        xPattDef = cosd(alphPatt).*gDef*(distSpaceT);
        yPattDef = sind(alphPatt).*gDef*(distSpaceT);
        % display ARP on the UE position probability map
        figure(figNumber); 
        subplot(1,2,subFigNumber); 
        [X,Y] = meshgrid(-25:25, 0:Due+15);
        p = mvnpdf([X(:) Y(:)], ueCoord([2,1]),...
            diag([stdCoords,stdCoords].^2));
        p = reshape(p,size(X));
        pcolor(X,Y,p); shading interp; hold on;
        plot(xPatt, yPatt, 'Color', '#D95319', 'LineWidth', 1.3);
        plot(xPattDef, yPattDef, 'Color', '#7E2F8E', 'LineWidth', 1.3);
        c = colorbar; c.Label.String = 'p(x,y)';
        axis tight;
        text(ueCoord(2), ueCoord(1), 'UE', 'HorizontalAlignment',...
            'center', 'VerticalAlignment', 'bottom')
        grid on; axis equal; xlabel('x, m'); ylabel('y, m');
        if ww == 0 
            legend('Gauss ARP','Location', 'southeast');
        elseif ww == 2
            legend('Rectangular ARP','Location', 'southeast');
        end
        subFigNumber = subFigNumber + 1;
        ueRxPwrPlt = [ueRxPwrPlt; ueRxPwr];
    end
    sgtitle([sprintf('RMSE %i m', s), antTypeCmt(antType)]);
    figNumber = figNumber + 1;
end
figure; plot(ueRxPwrPlt([1, 4, 2],:).')
grid on; xlabel('Calculation point number'); ylabel('P, dB');
legend('Without beam control', 'Rectangular ARP', 'Gauss ARP',...
    'Location', 'southeast');