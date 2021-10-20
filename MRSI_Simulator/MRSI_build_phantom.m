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
%   MemoryOptions: Key value argument pairs
%       Key: 'use_disc' Value: logical (1 or 0) 
%       Description: Boolean value to decide weather to save phantom object to disc as a binary file.
%       Use this option if your computer does not have enough memory.
%

%Output
%   phantom: a 2d matrix representing a phantom. phantom(i,j) represents a voxel.
%       phantom.x:     x coordinate in mm
%       phantom.y:     y coordinate in mm
%       phantom.met(:   vector of metabolites
%           phantom(i,j).met(j) :metabolite basis from sim_hamiltonian. Also includes density matrix

function [phantom] = MRSI_build_phantom(phan_size, met, b0, T2_star, MemoryOptions)
arguments
    phan_size (1,2) double = [200, 200]
    met (:, :) cell = make_cell()
    b0 (1,1) double {mustBeNonnegative} = 3
    T2_star (:,1) double = 0.1;
    MemoryOptions.use_disc (1,1) logical = 0
end
%initalize metabolites to be 46 long
counter = 1;


%get a list of all the metabolites that exist in the phantom
bool_metabolites = zeros(size(met,1), size(met,2), 46, 'logical');
for y = 1:size(met, 1) 
    for x = 1:size(met, 2)
        for m = 1:length(met{y,x})
            if(exist('all_mets', 'var'))
            %get list of names already found
                names = {all_mets.name};
            else
                names = '';
            end
            %if met name is not in list, add to list
            if(~strcmp(names, met{y,x}(m).name))
                all_mets(counter) = met{y,x}(m);

                all_mets(counter).shifts = all_mets(counter).shifts - 4.65;
                counter = counter + 1;
            end
            names = {all_mets.name};
            %logical array for position of metabolites
            bool_metabolites(y,x,strcmp(names, met{y,x}(m).name)) = 1;
        end
    end
end

%remove empty spaces in array
for m = 1:length(all_mets)
    if(strcmp(all_mets(m).name, ''))
        all_mets(m:end) = [];
        bool_metabolites(:,:,m:end) = [];
        break;
    end
end

if(length(T2_star) == 1)
    T2_star = repmat(T2_star, [length(all_mets), 1]);
end
%loop through all metabolites and get hamiltonian and density matrix
for m = length(all_mets):-1:1
    [phantom.met(m), d] = sim_Hamiltonian(all_mets(m), b0);
    phantom.d{m} = single(d{1});
    phantom.met_names{m} = all_mets(m).name;
    phantom.spins(m) = {zeros(size(met, 1), size(met, 2), size(d{1}, 1), size(d{1}, 2), 'single')};
    phantom.T2(m) = T2_star(m);
end
phantom.conc = zeros([size(met, [1,2]), length(all_mets)]);
%initalize phantom with empty structs
for m = 1:length(all_mets)
    for y = 1:size(bool_metabolites,1)
        for x = 1:size(bool_metabolites,2)
            if(bool_metabolites(y,x,m))
                
                %find all names of metabolites at that voxel
                names = {met{y,x}.name};
                %get index of metabolite currently looping over
                idx = strcmp(names, phantom.met_names{m});
                cur_met = met{y,x}(idx);
                if(isfield(cur_met, 'conc'))
                    
                    %multiply density matrix by concentration 
                    phantom.spins{m}(y,x,:,:) = phantom.d{m}*met{y,x}(idx).conc;
                    phantom.conc(y,x,m) = met{y,x}(idx).conc;
                else
                    %no concentration, assume uniform concentration
                    phantom.spins{m}(y,x,:,:) = phantom.d{m};
                    phantom.conc(y,x,m) = 1;
                end
            end
        end
    end
    if(MemoryOptions.use_disc)
        phantom.file{m} = save_spins(phantom.spins{m}, phantom.met_names{m}, 'MRSI_build_phantom');
        phantom.spins{m} = [];
    end
end

for m = 1:length(phantom.spins)
    phantom.spins{m} = complex(phantom.spins{m});
end


%set phantom phan_dY = phan_size(1)/size(met, 1);
phan_dX = phan_size(2)/size(met, 2);
phan_dY = phan_size(1)/size(met, 1);
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

%add x and y values to phantom
phantom.x = phan_x;
phantom.y = phan_y;

%if mets are empty fill with template values
if(isempty(all_mets))
    phantom.met = [];
    phantom.d ={0};
    phantom.met_names = {''};
    phantom.spins = {zeros(length(phantom.y), length(phantom.x))};
    phantom.T2 = T2_star;
end

end
function metabolites = make_cell()
load H2O.mat sysH2O;
metabolites = cell(32,32);
metabolites(10:12, 10:12) = {sysH2O};
end
