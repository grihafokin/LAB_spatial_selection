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
        stdCoordsArr = s; % RMSE of UE coordinate estimates along, m
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
        prbA = [0.1, 0.5, 0.9];
        [~, ~, ppR] = get_prob(ueCoord(2), ueCoord(1), stdCoords, prbA);
        
        ylimnin=-2*stdCoords; ylimnax=ppR(end) + ueCoord(1) + abs(ylimnin);
        xlimmin=5*stdCoords;
        [X,Y] = meshgrid(-xlimmin:xlimmin, ylimnin:ylimnax);
        p = mvnpdf([X(:) Y(:)], ueCoord([2,1]),...
            diag([stdCoords,stdCoords].^2));
        p = reshape(p,size(X)); pcolor(X,Y,p); hold on; shading interp;
        % c = colorbar; c.Label.String = 'p(x,y)'; 
        angA = 0:360;
        for i=1:length(prbA)
            plot(ppR(i)*cosd(angA) + ueCoord(2), ...
                ppR(i)*sind(angA) + ueCoord(1), 'k-');
            text(ueCoord(2), ppR(i) + ueCoord(1), num2str(prbA(i)), ...
           'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
        plot(xPatt,yPatt,'r--','LineWidth', 2.0); hold on;
        plot(xPattDef, yPattDef, 'g-', 'LineWidth', 1.5);
        text(ueCoord(2), ueCoord(1), 'UE', 'HorizontalAlignment',...
            'center', 'VerticalAlignment', 'bottom')
        grid on; xlabel('x, m'); ylabel('y, m'); 
        axis equal; axis([-xlimmin xlimmin ylimnin ylimnax-2]);
        if ww == 0 
            legend('Gaussian ARP','Location', 'southeast');
        elseif ww == 2
            legend('Rectangular ARP','Location', 'southeast');
        end
        subFigNumber = subFigNumber + 1;
        ueRxPwrPlt = [ueRxPwrPlt; ueRxPwr];
    end
    sgtitle([strcat(antTypeCmt(antType), sprintf('; RMSE = %i m',s))]);
    figNumber = figNumber + 1;
end
figure; 
plot(ueRxPwrPlt(1,:).', '-'); hold on;
plot(ueRxPwrPlt(2,:).', 'x-');
plot(ueRxPwrPlt(4,:).', 'o-');
% plot(ueRxPwrPlt([1, 4, 2],:).');
grid on; xlabel('Calculation point number'); ylabel('P, dB');
sgtitle([strcat(antTypeCmt(antType), sprintf('; RMSE = %i m',s))]);
legend('Without beam control', 'Rectangular ARP', 'Gaussian ARP',...
    'Location', 'southeast');