% plot SIR & delta SIR dependence on the number of antenna array elements
% and RMSE for scenarios: 1) maximum and maximum and null ARP beam shape
% control; 2) Gauss beam width control; 3) Rectangular  beam width control;
% delta SIR dependence for all scenarios is calculated w.r.t. only maximum 
% ARP beam shape control (maximum ARP beam shape control is plotted dotted)
clear all; close all; clc;
stdCoordsArr = [1, 3, 10];  % RMSE of UE coordinate estimates, m
NellArr = 4:4:20;           % number of AA elements in one dimension
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
antPattCntrlArr = [0, 1, 3, 3];
win_typeArr = [0, 0, 0, 2];
antPattCntrlCmt = ["Maximum ARP beam shape control", ...
                   "Maximum & null ARP beam shape control", ...
                   "Gauss beam width control", ...
                   "Rectangular  beam width control"];
antTypeCmt = ["URA", "ULA", "UCA"];
Zaa = [];
for aa=1:length(antPattCntrlArr)
    Zst = [];
    antPattCntrl = antPattCntrlArr(aa);
    win_type = win_typeArr(aa);
    for s=stdCoordsArr
        Zn = [];
        stdCoords = s;
        for nel=NellArr
            Zsc = [];
            Nel = nel;
            switch win_type
                case 0
                    sc = 3.79*exp(-0.126*stdCoords);
                case 1
                    sc = 0.1;
                case 2
                    sc = -0.2*stdCoords + 4.65;
            end
            rng('default')
            fprintf('alg %i rmse=%5.1f, nel=%3i, scl=%4.2f\n',...
                antPattCntrl, s, nel, sc);
            [X,Y,Z] = lab_beam_shape_fcn(antType,Nel,stdCoords,...
                                            antPattCntrl,win_type);
            Zsc = [Zsc, mean(Z(:))];
            Zn = [Zn;Zsc];
        end
        Zst{end+1} = Zn;
    end
    Zaa{end+1} = Zst;
end
% plotting figures
ZZmS_1 = Zaa{1};
for aa=2:length(antPattCntrlArr)
    ZZmS = Zaa{aa};
    figure; hold on; grid on; axis tight;
    for kk=1:size(ZZmS,2)
        plot(NellArr, ZZmS_1{kk}, '-');
    end
    set(gca,'ColorOrderIndex',1)
    for kk=1:size(ZZmS,2)
        plot(NellArr, ZZmS{kk}, '--');
    end
    xlabel('N'); ylabel('SIR_{avg}, dB');
    legFrst = ['RMSE = ',num2str(stdCoordsArr(1)), ' m'];
    if antPattCntrlArr(aa) == 3
        legend([legFrst, strcat(string(stdCoordsArr(2:end)),' m')],...
            'NumColumns',3, 'Location', 'southeast');
    else
        legend([legFrst, strcat(string(stdCoordsArr(2:end)),' m')],...
            'NumColumns',3);
    end
    title([antPattCntrlCmt(aa), antTypeCmt(antType)]);
    
    figure; hold on; grid on; axis tight;
    for kk=1:size(ZZmS,2)
        plot(NellArr, ZZmS{kk}-ZZmS_1{kk}, '-');
    end
    xlabel('N'); ylabel('\DeltaSIR_{avg}, dB');
    legFrst = ['RMSE = ',num2str(stdCoordsArr(1)), ' m'];
    if antPattCntrlArr(aa) == 3
        legend([legFrst, strcat(string(stdCoordsArr(2:end)),' m')],...
            'NumColumns',3, 'Location', 'southeast');
    else
        legend([legFrst, strcat(string(stdCoordsArr(2:end)),' m')],...
            'NumColumns',3);
    end
    title([antPattCntrlCmt(aa), antTypeCmt(antType)]);
end