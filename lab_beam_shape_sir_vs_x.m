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
Nel = 8;  
% RMSE of UE coordinate estimates along the x, y, z axes in meters
stdCoordsArr = [0, 5]; 
% selection of the NB beamforming (BF) control algorithm:
%   0 - maximum antenna radiation pattern (ARP) beam shape control
%   1 - maximum and null antenna radiation pattern (ARP) beam shape control
%   2 - adaptive antenna radiation pattern (ARP) beam shape control
%   3 - antenna radiation pattern (ARP) beam width (BW) control
antPattCntrlArr = [0, 1, 2];

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
    for ii=1:2
        plot(X, ZZ(ii,:),'--');
    end
    grid on; xlabel('x, m'); ylabel('SIR, dB');
    legend('maximum beam shape control',...
           'maximum and null beam shape control');
end
figure; hold on
for ii=1:2
    plot(X, ZZa(ii,:),'--');
end
grid on; xlabel('x, m'); ylabel('SIR, dB');
legend(sprintf("RMSE = %i m",stdCoordsArr(1)),...
       sprintf("RMSE = %i m",stdCoordsArr(2)));