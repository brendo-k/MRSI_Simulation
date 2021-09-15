%MRSI_build_phantom.m
%This function creates a phantom that can be then used for simulations. THe
%scales are inverted because matlab uses the first dimension as y and the
%second as x (This can be seen if you imagine matrix indexing).
%
%
%Inputs
%   phan_size:  a 2 element vector of the size in y and x direction in mm
%               (eg. 200, 200 is 200 mm in both y and x directions)
%   met:        2d cell array of metabolites. met(i,j) represents metabolites at
%               position y(i), x(j)
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
%initalize metabolites to be 46 long
all_mets(1:46) = struct('J', 0, 'shifts', 0, 'name', '', 'scaleFactor', 0);
counter = 1;
bool_metabolites = zeros(size(met,1), size(met,2), 46);
for y = 1:size(met, 1) 
    for x = 1:size(met, 2)
        for m = 1:length(met{y,x})
            names = {all_mets.name};
            if(~strcmp(names, met{y,x}(m).name))
                all_mets(counter) = met{y,x}(m);
                all_mets(counter).shifts = all_mets(counter).shifts - 4.65;
                counter = counter + 1;
            end
            bool_metabolites(y,x,strcmp(names, met{y,x}(m).name)) = 1;
        end
    end
end
bool_metabolites = logical(bool_metabolites);
for m = 1:length(all_mets)
    if(strcmp(all_mets(m).name, ''))
        all_mets(m:end) = [];
        break;
    end
end

for m = length(all_mets):-1:1
    [phantom.met(m), d] = sim_Hamiltonian(all_mets(m), b0);
    phantom.d{m} = d{1};
    phantom.met_names{m} = all_mets(m).name;
    phantom.spins(m) = {zeros(size(met, 1), size(met, 2), size(d{1}, 1), size(d{1}, 2), 'single')};
end

names = {all_mets.name};
%initalize phantom with empty structs
for m = 1:length(all_mets)
    for y = 1:size(bool_metabolites,1)
        for x = 1:size(bool_metabolites,2)
            phantom.spins{m}(y,x,:,:) = phantom.d{m};
        end
    end
end


%set phantom parameters
phan_dY = phan_size(1)/size(met, 1);
phan_dX = phan_size(2)/size(met, 2);

%x coordinates for the phantom
phan_y = -phan_size(1)/2 + phan_dY/2:phan_dY:phan_size(1)/2 - phan_dY/2;
phan_x = -phan_size(2)/2 + phan_dX/2:phan_dX:phan_size(2)/2 - phan_dX/2;

%if the y dimension is empty allow one corrdinate of zero
if(isempty(phan_y))
    phan_y = 0;
end

%if the x dimension is empty allow one corrdinate of zero
if(isempty(phan_x))
    phan_x = 0;
end

phantom.x = phan_x;
phantom.y = phan_y;
if(isempty(all_mets))
    phantom.met = [];
    phantom.d ={0};
    phantom.met_names = {''};
    phantom.spins = {zeros(length(phantom.y), length(phantom.x))};
end

end
function metabolites = make_cell()
load H2O.mat sysH2O;
metabolites = cell(128,128);
metabolites(10:12, 10:12) = {sysH2O};
end
