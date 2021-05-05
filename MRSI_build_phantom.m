%Inputs
%   nx: number of points in the x direction 
%   ny: number of points in the y direction
%   sizeX: the size of the phantom in the x direction
%   sizeY: the size of the phantom in the y direciton

%Output
%   phantom: a [nx, ny] array containing voxel object
%       phantom(i,j).d: density matrix at position i,j
%       phantom(i,j).dI: chemical shift at position i,j

function [phantom] = MRSI_build_phantom(nx, ny, sizeX, sizeY)
%set phantom parameters
phantom_deltaX = sizeX/nx;
phantom_deltaY = sizeY/ny;

%Set some bases TODO: get different bases for different metabolites
I0=complex([1 0;0 1]);
Ix=complex(0.5*[0 1;1 0]);
Iy=(1i/2)*[0 -1;1 0];
Iz=complex((1/2)*[1 0;0 -1]);
dI = 0;
d0 = Iz;

%x coordinates for the phantom
phan_x = -sizeX/2 + phantom_deltaX/2:phantom_deltaX:sizeX/2 - phantom_deltaX/2;
phan_y = -sizeY/2 + phantom_deltaY/2:phantom_deltaY:sizeY/2 - phantom_deltaY/2;
if(isempty(phan_y))
    phan_y = [0];
end
if(isempty(phan_x))
    phan_x = [0];
end

%initalize phantom with empty structs
phantom = struct('d', cell(length(phan_x), length(phan_y)), 'dI', cell(length(phan_x),length(phan_y)));
%create phantom
len_x = length(phan_x);
len_y = length(phan_y);
sig_x = ceil(len_x/4);
sig_y = ceil(len_y/4);

for i = 1:length(phan_x)
    for j = 1:length(phan_y)
        %Set coordinate
        phantom(i,j).x = phan_x(i);
        phantom(i,j).y = phan_y(j);
        
        %water from 20 to 44, already excited
        if i >= sig_x && i <= len_x - sig_x + 1 && j >= sig_y/2 && j <= len_y - sig_y/2 + 1
            phantom(i,j).d = d0;
            phantom(i,j).dI = dI;
%         elseif i > 40 && i <= 50
%             phantom(i,j).d = d0;
%             phantom(i,j).dI = dI + 1;
        else
            %Nothing elsewhere
            phantom(i,j).d = I0;
            phantom(i,j).dI = 0;
        end
    end
end
end