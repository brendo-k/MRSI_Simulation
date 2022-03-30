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

function [phantom] = MRSI_build_phantom(phantomSize, metabolites, concMap, b0, ...
                                        T2, MemoryOptions)
    arguments
        phantomSize (1, 2) double = [200, 200]
        metabolites (1, :) struct = make_cell()
        concMap (:, :, :) = ones(length(metabolites), 20, 20)
        b0 (1, 1) double {mustBeNonnegative} = 3
        T2 (1, :) double = 0.1;
        MemoryOptions.use_disc (1,1) logical = 0
    end

    if(length(metabolites) ~= size(concMap, 1))
        numMetabolites = getNumMetabolites(metabolites);
        if(numMetabolites ~= size(concMap, 1))
            error('Not enough metabolites in concentration Map')
        end
        
        newConcentrationMap = expandConcentrationMap(metabolites, concMap);
        concMap = newConcentrationMap;
    end

    [phantom, numXPoints, numYPoints] = calculateFovandVoxelSize(phantomSize, concMap);
    
    if(length(T2) == 1)
        T2 = repmat(T2, [1, length(metabolites)]);
    end

    %loop through all metabolites and get hamiltonian and density matrix
    for m = length(metabolites):-1:1
        [phantom.met(m), d] = sim_Hamiltonian(metabolites(m), b0);
        phantom.d{m} = complex(single(d{1}));
        phantom.met_names{m} = metabolites(m).name;
        phantom.spins{m} = zeros([size(d{1}), size(concMap, [2, 3])], 'single');
        phantom.T2(m) = T2(m);
    end

    phantom.conc = zeros([size(metabolites, [1, 2]), length(metabolites)]);
    %initalize phantom with empty structs
    for m = 1:length(metabolites)
        for y = 1:numYPoints
            for x = 1:numXPoints
                %multiply density matrix by concentration
                phantom.spins{m}(:, :, y, x) = phantom.d{m}*concMap(m, y, x);
            end
        end
        if(MemoryOptions.use_disc)
            phantom.file{m} = save_spins(phantom.spins{m}, phantom.met_names{m}, 'MRSI_build_phantom');
            phantom.spins{m} = [];
        end
    end
    phantom.conc = concMap;
end



function metabolites = make_cell()
    load H2O.mat sysH2O;

    metabolites = [sysH2O];
end

function [phantom, numX, numY] = calculateFovandVoxelSize(phantomSize, concMap, phantom)
    %Creates cartesian grid of coordinates over the phantom
    phantomSizeY = phantomSize(1); phantomSizeX = phantomSize(2);
    [numY, numX] = size(concMap, [2,3]);
    %set voxel size in teh x and y direction
    deltaX = phantomSizeX/numX;
    deltaY = phantomSizeY/numY;
    % calculate coordinates for the phantom
    yCoordinates = createCoordinates(phantomSizeY, deltaY);
    xCoordinates = createCoordinates(phantomSizeX, deltaX);

    %add x and y values to phantom
    phantom.x = xCoordinates;
    phantom.y = yCoordinates;
    phantom.voxelSizeX = deltaX;
    phantom.voxelSizeY = deltaY;
end

function numMetabolites = getNumMetabolites(metabolites)
    metNames = {};
    for iMet = 1:length(metabolites)
        metaboliteName = exractMetaboliteName(metabolites(iMet).name);
        if(contains(metNames, metaboliteName))
            metNames{end + 1} = metaboliteName;
        end
    end
    numMetabolites = length(metNames);
end

function newConcentrationMap = expandConcentrationMap(metabolites, concMap)
    newConcentrationMapSize = [length(metabolites), size(concMap, [2, 3])];
    newConcentrationMap = zeros(newConcentrationMapSize);
    prevMetabolite = extractMetaboliteName(metabolites(1).name);
    metaboliteCounter = 1;
    for iMetabolite = 1:length(metabolites)
        currentMetabolite = extractMetaboliteName(metabolites(iMetabolite).name);
        if(~isequal(currentMetabolite, prevMetabolite))
            metaboliteCounter = metaboliteCounter + 1;
            newConcentrationMap(iMetabolite, :, :) = concMap(metaboliteCounter, :, :);
        else
            newConcentrationMap(iMetabolite, :, :) = concMap(metaboliteCounter, :, :);
        end
    end
end