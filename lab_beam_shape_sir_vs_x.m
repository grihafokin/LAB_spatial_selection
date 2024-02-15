% plot SIR dependence on distance for scenario: 1) comparison of maximum 
% and  maximum and null ARP beam shape control for RMSE=0; 2) comparison 
% of maximum and  maximum and null ARP beam shape control for RMSE=5; 
% 3) adaptive ARP beam shape control
clear all; close all; clc; 
% selection of antenna array type
% 1 - planar or uniform rectangural antenna array (URA)
% 2 - uniform linear antenna array (ULA)
% 3 - uniform circular antenna array (UCA); supports BF: antPattCntrl=0,1,2
antType = 1;
% number of AA elements in one dimension
Nel = 20;  
% RMSE of UE coordinate estimates along the x, y, z axes in meters
stdCoordsArr = [0 10]; 
% selection of the NB beamforming (BF) control algorithm:
%   0 - maximum antenna radiation pattern (ARP) beam shape control
%   1 - maximum and null antenna radiation pattern (ARP) beam shape control
%   2 - adaptive antenna radiation pattern (ARP) beam shape control
%   3 - antenna radiation pattern (ARP) beam width (BW) control
antPattCntrlArr = [0, 1];

ZZa = [];
for s=1:length(stdCoordsArr)
    ZZ = [];
    stdCoords = stdCoordsArr(s);
    for aa=antPattCntrlArr
        antPattCntrl = aa;
        if aa ~= 2
            rng('default');
        else
            rng(s);
        end
        [X,Y,Z] = lab_beam_shape_fcn(antType,Nel,stdCoords,antPattCntrl);
        ZZ = [ZZ; Z];
    end
    ZZa = [ZZa; Z];
    figure; hold on;
    for ii=1:length(antPattCntrlArr)
        plot(X, ZZ(ii,:),'-', 'linewidth',2);
    end
    grid on; xlabel('x, m'); ylabel('Instantaneous SIR, dB'); axis('tight');
    legend('maximum beam shape control',...
           'maximum and null beam shape control');
    title(sprintf("URA %i×%i, RMSE = %i m",Nel, Nel, stdCoords));
end
figure; hold on; axis('tight');
for ii=1:size(ZZa,1)
    plot(X, ZZa(ii,:),'-','linewidth',2);
end
grid on; xlabel('x, m'); ylabel('SIR, dB');
title(sprintf("Joint maximum and null beam shape control, URA %i×%i",Nel, Nel));
legend(sprintf("RMSE = %i m",stdCoordsArr(1)),...
       sprintf("RMSE = %i m",stdCoordsArr(2)));
axis('tight');