% 
clear all; close all; clc;
stdCoordsArr = [1, 3, 5, 10];  % RMSE of UE coordinate estimates, m
NellArr = 4:4:20;              % number of AA elements in one dimension
% selection of antenna array type
% 1 - planar or uniform rectangural antenna array (URA)
% 2 - uniform linear antenna array (ULA)
% 3 - uniform circular antenna array (UCA); supports BF: antPattCntrl=0,1,2
antType = 1;
% vector of the NB beamforming (BF) control algorithm:
%   0 - maximum antenna radiation pattern (ARP) beam shape control
%   1 - maximum and null antenna radiation pattern (ARP) beam shape control
%   2 - adaptive antenna radiation pattern (ARP) beam shape control
%   3 - antenna radiation pattern (ARP) beam width (BW) control; adaptive 
%     beam shape control is out of scope, because does not depend on RMSE
antPattCntrlArr = [3, 3];
% select ARP beam shape for beam width control algorithm (antPattCntrl = 3)
% 0 - Gauss window
% 1 - raised cosine window
% 2 - rectangular window
win_typeArr = [2, 0];
antTypeCmt = ["URA", "ULA", "UCA"];
Zaa = [];
for aa=1:length(antPattCntrlArr)
    Zst = [];
    antPattCntrl = antPattCntrlArr(aa);
    win_type = win_typeArr(aa);
    switch win_type
        case 0
            sclArr = 0.1:0.5:5.5; % Gauss window
        case 1
            sclArr = 0:0.1:0.9;   % Raised cosine window; have not been sim
        case 2
            sclArr = 0.5:0.5:5.5; % Rectangular window
    end
    for s=stdCoordsArr
        Zn = [];
        stdCoords = s;
        for nel=NellArr
            Zsc = [];
            Nel = nel;
            for scl=sclArr
                rng('default')
                fprintf('win %i std=%.1f, Nel=%i, Scl=%.2f\n', ...
                    win_type, s, nel, scl);
                [X,Y,Z] = lab_beam_shape_fcn(antType,Nel,stdCoords,...
                            antPattCntrl,win_type,scl);
                Zsc = [Zsc, mean(Z(:))];
            end
            Zn = [Zn;Zsc];
        end
        Zst{end+1} = Zn;
    end
    Zaa{end+1} = Zst;
end

% plotting results
for aa=1:length(antPattCntrlArr)
    figure;
    ZZmS = Zaa{aa};
    [X,Y] = meshgrid(NellArr, sclArr);
    for kk=1:size(ZZmS,2)
        subplot(2,2,kk)
        surf(X, Y, ZZmS{kk}.', 'FaceColor', 'interp', 'EdgeColor','none');
        c1 = colorbar; c1.Label.String = 'SIR_{avg}, dB';
    grid on; xlabel('N'); ylabel('s'); view([0, 90]); axis tight
    legend(sprintf('RMSE = %d m',stdCoordsArr(kk)),'Location','southeast');
    title(sprintf('RMSE = %d m',stdCoordsArr(kk)));
    end
    if win_typeArr(aa) == 0
        sgtitle(sprintf('Gauss window ARP %s', antTypeCmt(antType)));
    elseif win_typeArr(aa) == 2
        sgtitle(sprintf('Rectangular window ARP %s', antTypeCmt(antType)));
    end
end