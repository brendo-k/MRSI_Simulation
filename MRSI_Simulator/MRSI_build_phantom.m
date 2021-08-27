%MRSI_build_phantom.m
%This function creates a phantom that can be then used for simulations
%
%
%Inputs
%   phan_size:  a 2 element vector of the size in x and y direction in mm
%               (eg. 200, 200 is 200 mm in both x and y directions)
%   met:        2d cell array of metabolites. met(i,j) represents metabolites at
%               position x(i), y(j)
%   b0:         b0 field strength
%

%Output
%   phantom: a 2d matrix representing a phantom. phantom(i,j) represents a voxel.
%       phantom(i,j).x:     x coordinate in mm
%       phantom(i,j).y:     y coordinate in mm
%       phantom(i,j).met:   vector of metabolites
%           phantom(i,j).met(j) :metabolite basis from sim_hamiltonian. Also includes density matrix

function [phantom] = MRSI_build_phantom(phan_size, met, b0)
arguments
    phan_size (1,2) double = [200, 200]
    met (:, :) cell = make_cell()
    b0 = 3;
end

%set phantom parameters
phan_dX = phan_size(1)/size(met, 1);
phan_dY = phan_size(2)/size(met, 2);

%x coordinates for the phantom
phan_x = -phan_size(1)/2 + phan_dX/2:phan_dX:phan_size(1)/2 - phan_dX/2;
phan_y = -phan_size(2)/2 + phan_dY/2:phan_dY:phan_size(2)/2 - phan_dY/2;


%make sure the size of phantom matches the size of met cell array
assert(length(phan_x) == size(met, 1), 'Length of phantom along x should equal length of met along first dimension')
assert(length(phan_y) == size(met, 2), 'Length of phantom along y should equal length of met along first dimension')

%if the y dimension is empty allow one corrdinate of zero
if(isempty(phan_y))
    phan_y = 0;
end

%if the x dimension is empty allow one corrdinate of zero
if(isempty(phan_x))
    phan_x = 0;
end

%initalize phantom with empty structs
for j = length(phan_x):-1:1
    for i = length(phan_y):-1:1
        %Set coordinate
        phantom(i,j).x = phan_x(j);
        phantom(i,j).y = phan_y(i);
        if(isempty(met{i,j}))
            phantom(i,j).met = [];
            phantom(i,j).d = [];
            continue
        end
        for m = 1:length(met{i,j})
            met{i,j}(m).shifts = met{i,j}(m).shifts - 4.65;
        end
        [metabolites, d] = sim_Hamiltonian(met{i,j}, b0);
        phantom(i,j).met = metabolites;
        phantom(i,j).d = d;
    end
end
end

function metabolites = make_cell()
load H2O.mat sysH2O;
metabolites = cell(16,16);
metabolites(10:12, 10:12) = {sysH2O};
end
