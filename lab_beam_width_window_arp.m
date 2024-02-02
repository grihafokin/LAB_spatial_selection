clear all; close all; clc;
c = physconst('LightSpeed');
f = 30e9;       % carrier frequency in mmWave (millimetre wave) range, Hz
lamb = c/f;     % wavelength, m
da = 0.5*c/f;   % distance between antenna array (AA) elements, m
BW = 50;        % main lobe ARP beam width, degrees
Nel = 20;       % number of AA elements in one dimension
% selection of antenna array type
% 1 - planar or uniform rectangural antenna array (URA)
% 2 - uniform linear antenna array (ULA)
% 3 - uniform circular antenna array (UCA); is not supported
antType = 1;
antElPos = createAnt(antType, Nel, da); % construction of AA

figNumber = 1;
for j=[0,1]
    % calculation of coefficient for rectangular window ARP
    [wr, azAngP, antPattPr] = beamshapingWeight(2, BW, 0, Nel, 1, j);
    if (antType == 1)
        % calculation of the vector of weighting 
        % coefficients of URA for vertical AA elements
        wr = repmat(wr, Nel, 1)/Nel;
        wr = wr(:);
    end
    gr = zeros(1, length(azAngP));
    for i=1:length(azAngP)
        gr(i) = getAntPatternG(antElPos, f, azAngP(i), 0, wr, 0);
    end

    figure;
    plot(azAngP, antPattPr,'LineWidth', 2); hold on; 
    plot(azAngP, gr,'--','LineWidth', 2); grid on;
    xlabel('\phi, \circ'); ylabel('|A(\phi)|');
    legend('given window', 'synthesized');
    figNumber = figNumber + 1;
end
% calculation of coefficient for Gauss window ARP
[wg, azAngP, antPattPg] = beamshapingWeight(0, BW, 0, Nel, 1);
if (antType == 1)
    % calculation of the vector of weighting 
    % coefficients of URA for vertical AA elements
    wg = repmat(wg, Nel, 1)/Nel;
    wg = wg(:);
end
gg = zeros(1, length(azAngP));
for i=1:length(azAngP)
    gg(i) = getAntPatternG(antElPos, f, azAngP(i), 0, wg, 0);
end

figure;
plot(azAngP, antPattPg, 'LineWidth', 2); hold on;
plot(azAngP, gg, '--', 'LineWidth', 2); grid on;
xlabel('\phi, \circ'); ylabel('|A(\phi)|');
legend('given window', 'synthesized');