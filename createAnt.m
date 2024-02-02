% antenna array creation function 
% input parameters:
%   antType - antenna arrayА type
%   Nel     - number of antenna array elements
%   d       - distance between antenna array elements, m
% output parameters:
%   antPos  - antenna array elements coordinates [x,y,z], m
%            (size [Nel x 3] for ULA, [Nel^2 x 3] for URA)
function antPos = createAnt(antType, Nel, d)
switch antType
    case 1 % Planar or Uniform Rectangular Array (URA) 
        % number of antenna array elements on one side of the square
        % of antenna array; total number of elements Nel^2
        yCoords = repmat((-(Nel-1)/2:1:(Nel-1)/2).'*d, Nel, 1);
        zCoords = reshape(repmat((-(Nel-1)/2:1:(Nel-1)/2)*d, Nel, 1),[],1);
        zCoords = zCoords(:);
        xCoords = zeros(size(yCoords));
        antPos = [xCoords, yCoords, zCoords];        
    case 2 % Uniform Linear Array (ULA) 
        % Uniform Linear antenna array made of Nel elements
        antLen = (Nel-1)*d;
        yCoords = (-antLen/2:d:antLen/2).';
        zCoords = zeros(size(yCoords));
        xCoords = zeros(size(yCoords));
        antPos = [xCoords, yCoords, zCoords];  
    case 3 % Uniform Circular Array (UCA)  
        % Uniform Circular antenna array made of Nel = Nel*2 elements
        % calculation of the radius of circular lattice using the arc 
        % chord length formula; the chord length is equal to the specified 
        % distance between the antenna array elements d
        r = d/2/sind(360/Nel/2);
        antDph = (0:360/Nel:360 - 360/Nel).';
        xCoords = r*cosd(antDph);
        yCoords = r*sind(antDph);
        zCoords = zeros(size(yCoords));
        antPos = [xCoords, yCoords, zCoords];  
end
end