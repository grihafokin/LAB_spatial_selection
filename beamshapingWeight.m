% function for calculating the vector of weighting coefficients for antenna
% radiation pattern (ARP) synthesis of given main beam lobe shape and width

% input parameters:
%   win_type - type of window function for synthesizing ARP
%              (rectangular, raised cosine, Gaussian)
%   BW -       ARP main beam lobe width, degrees
%   stAng -    direction of the ARP maximum in the azimuthal plane, degrees
%   scl -      ARP main beam lobe width multiplier
% output parameters:
%   w -        vector of weight coefficients
%   antPattP - given form of ARP
%   azAngP -   values of azimuthal angles in which antPattP is set
function [w, azAngP, antPattP]  = ...
    beamshapingWeight(win_type, BW, stAng, N, scl, varargin)
k = 0:N-1;
% calculation of the points at which the synthesized ARP will be specified
psi_k = wrapToPi((k - (N - 1)/2)*2*pi/N);
% conversion of psi_k to azimuth angles; for d=lambda/2 formula simplified 
% 2*pi/lamb*sin(az)*d = pi*sin(az); sign "-" for ease of data presentation
azAngW = -asind(psi_k/pi);
Na = length(azAngW);
azAngP = -90:0.1:89;
Np = length(azAngP);
antPattP = zeros(1, Np);
% formation of ARP of a given shape and width
switch win_type
    case 0  % ARP in the form of a Gaussian function
        antPattW = exp(-((azAngW-stAng).^2/2/(BW*scl/2)^2));
        antPattP = exp(-((azAngP-stAng).^2/2/(BW*scl/2)^2));
    case 1 % ARP in form of raised cosine function; not used in simulation
        beta = 0.35;
        B = BW*scl/(1-beta);
        antPattW = zeros(1, Na);
        for j=1:Na
            if (azAngW(j) > -(1-beta)*B/2+stAng && ...
                    azAngW(j) < (1-beta)*B/2+stAng)
                antPattW(j) = 1;
            elseif (azAngW(j) > (1-beta)*B/2+stAng &&...
                    azAngW(j) < (1+beta)*B/2+stAng) ||...
                    (azAngW(j) < -(1-beta)*B/2+stAng &&...
                    azAngW(j) > -(1+beta)*B/2+stAng)
                antPattW(j) = 0.5*(1+cos(pi/beta/B*(abs(azAngW(j)-stAng)...
                    - (1-beta)*B/2)));
            end
        end
        for j=1:Np
            if (azAngP(j) > -(1-beta)*B/2+stAng &&...
                    azAngP(j) < (1-beta)*B/2+stAng)
                antPattP(j) = 1;
            elseif (azAngP(j) > (1-beta)*B/2+stAng &&...
                    azAngP(j) < (1+beta)*B/2+stAng) ||...
                    (azAngP(j) < -(1-beta)*B/2+stAng &&...
                    azAngP(j) > -(1+beta)*B/2+stAng)
                antPattP(j) = 0.5*(1+cos(pi/beta/B*(abs(azAngP(j)-stAng)...
                    - (1-beta)*B/2)));
            end
        end
    case 2 % rectangular ARP
        antPattW = zeros(size(azAngW));
        antPattW(azAngW>(-BW*scl/2+stAng) & azAngW<(BW*scl/2+stAng)) = 1;
        % additional handling of the case, when there is only one nonzero 
        % value in the antPattW; in this case, the neighboring points on 
        % the left and on the right are interpolated linearly
        rEdge = (-BW*scl/2+stAng);
        lEdge = (BW*scl/2+stAng);
        if (sum(antPattW) == 1)
            lInd = find(azAngW>lEdge,1,'last');
            rInd = find(azAngW<rEdge,1,'first');
            if (lInd == 1)
                antPattW(lInd) = 1 - (azAngW(lInd) - lEdge)/...
                    (azAngW(lInd) - azAngW(lInd+1))/2;
            elseif (lInd > 1)
                antPattW(lInd) = 1 - (azAngW(lInd) - lEdge)/...
                    (azAngW(lInd-1) - lEdge);
            end
            if (rInd == N)
                antPattW(rInd) = 1 - (azAngW(rInd) - rEdge)/...
                    (azAngW(rInd) - azAngW(rInd-1))/2;
            elseif (rInd < N)
                antPattW(rInd) = 1 - (azAngW(rInd) - rEdge)/...
                    (azAngW(rInd+1) - rEdge);
            end
        end
        antPattP(azAngP>(-BW*scl/2+stAng) & azAngP<(BW*scl/2+stAng)) = 1;
end 
Ak = conj(antPattW).*exp(-1i*psi_k*(N-1)/2); % calculation of Ak
n = 0:N-1;
% calculation of bn; ifft is used to speed up the calculation
bn = ifft(Ak); % bn = Ak*exp(1i*(k'*n)*2*pi/N)/N;
% calculation of vector of coefficients
w = (bn.*exp(-1i*n*pi*(N-1)/N))';
% additional processing of coefficient vector w by Hamming window function
if (nargin == 6)
    useW = varargin{1};
else
    useW = 1;
end
if win_type == 2 && useW == 1
    alph = 0.54;
    % Hamming window coefficients
    wW = alph - (1-alph)*cos(2*pi*(0:N-1)/(N-1));
    w = w.*wW.';
end
end