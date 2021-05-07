%Inputs
%   n: number of points in the x and y direction
%   size: Size in x and y direction
%   lb: lower bounds of the phantom in x and y coordinates

%Output
%   phantom: a [nx, ny] array containing voxel object
%       phantom(i,j).d: density matrix at position i,j
%       phantom(i,j).dI: chemical shift at position i,j

function [phantom] = MRSI_build_phantom(n, size, lb, ub, met)
%tests to make sure input arguments are vectors of size 1,2
validateattributes(n,{'numeric'},{'size', [1,2]})
validateattributes(size,{'numeric'},{'size', [1,2]})
validateattributes(lb,{'numeric'},{'size', [1,2]})
validateattributes(ub,{'numeric'},{'size', [1,2]})

%set phantom parameters
phan_dX = size(1)/n(1);
phan_dY = size(2)/n(2);

%Set some bases TODO: get different bases for different metabolites
I0=complex([1 0;0 1]);
Iz=complex((1/2)*[1 0;0 -1]);
dI = 1;
d0 = Iz;

%x coordinates for the phantom
phan_x = -size(1)/2 + phan_dX/2:phan_dX:size(1)/2 - phan_dX/2;
phan_y = -size(2)/2 + phan_dY/2:phan_dY:size(2)/2 - phan_dY/2;

if(isempty(phan_y))
    phan_y = 0;
end
if(isempty(phan_x))
    phan_x = 0;
end

%initalize phantom with empty structs
phantom = struct('d', cell(length(phan_x), length(phan_y)), 'dI', cell(length(phan_x),length(phan_y)));

for i = 1:length(phan_x)
    for j = 1:length(phan_y)
        %Set coordinate
        phantom(i,j).x = phan_x(i);
        phantom(i,j).y = phan_y(j);
        
        %metabolites from lb to ub
        if i >= lb(1) && i <= ub(1) && j >= lb(2) && j <= ub(2)
             phantom(i,j).d = d0;
             phantom(i,j).dI = dI;
        else
            %Nothing elsewhere
            phantom(i,j).d = I0;
            phantom(i,j).dI = 0;
        end
    end
end
end