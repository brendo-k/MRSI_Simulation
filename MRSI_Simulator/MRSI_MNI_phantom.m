function [phantom] = MRSI_MNI_phantom(B0, slice, mets, resolution, MemoryOptions)
arguments
    B0 (1, 1) double
    slice (1, 1) double {mustBePositive}
    mets (:, 1) cell = [];
    resolution.resolution (1, :) {mustBeMember(resolution.resolution, {'1mm', '2mm'})} = '2mm'
    MemoryOptions.use_disc (1,1) logical = 0
end
    [whiteMatterIntensity, greyMatterIntensity, protonDensity, skullIntensity, ...
        whiteMatterCoordinates] = readAtlas(slice, resolution.resolution);
    
    mniAtlasSize = size(whiteMatterIntensity);
    phantom = calculateSpatialParameters(mniAtlasSize, whiteMatterCoordinates);
    
    % start adding metabolites
    % if input metabolites are empty, add all metabolites
    if(isempty(mets))
        %load spin system
        spins = load('spinSystems.mat');
        lipids = load('Lip.mat');
        metabolites = [spins.sysCr, spins.sysNAA, spins.sysNAAG,...
            spins.sysPCh, spins.sysGPC, spins.sysIns, spins.sysGlu,...
            spins.sysGln, spins.sysH2O, lipids.sysLip];
    else
        metabolites = mets;
    end
    
    % metabolite labels
    [whiteMatterConcDict, greyMatterConcDict] = getConcentrations;
    numMetabolites = length(metabolites);
    
    % Cell array of metabolite names that are loaded in
    metaboliteNames = {metabolites.name};
    % Fid-A loads some metabolites in multiple parts. Count the occurences of
    % the multiple parts and scale down the concentration accordingly
    metaboliteLabels = greyMatterConcDict.keys;
    for iKey = 1:length(metaboliteLabels)
        currentLabel = metaboliteLabels{iKey};
        if(strcmp(currentLabel, 'Cho'))
            numMetaboliteParts = sum(contains(metaboliteNames, 'PCh'));
            numMetaboliteParts = numMetaboliteParts + sum(contains(metaboliteNames, 'GPC'));
        else
            numMetaboliteParts = sum(contains(metaboliteNames, currentLabel));
        end
        
        %scale down the concentrations
        whiteMatterConcDict(currentLabel) = whiteMatterConcDict(currentLabel)/numMetaboliteParts;
        greyMatterConcDict(currentLabel) = greyMatterConcDict(currentLabel)/numMetaboliteParts;
    end
    

    concentrationMap = zeros([numMetabolites, mniAtlasSize]);
    for iMetabolite = 1:numMetabolites
        currentMetabolite = metabolites(iMetabolite);
        metaboliteLabel = currentMetabolite.name;
        
        metaboliteName = extractMetaboliteName(metaboliteLabel);
        
        if(contains(metaboliteName, 'lipids', 'IgnoreCase', true))
            concentrationMap(iMetabolite, :, :, :) = skullIntensity * 1000;
        elseif(contains(metaboliteName, 'h2o', 'IgnoreCase', true))
            concentrationMap(iMetabolite, :, :, :) = protonDensity;
        else
            concentrationMap(iMetabolite, :, :, :) = ... 
                 whiteMatterIntensity * whiteMatterConcDict(metaboliteName) + ...
                 greyMatterIntensity * greyMatterConcDict(metaboliteName);
        end
    end
    %build concentration map
    concentrationMap = concentrationMap(:, :, :, slice);

    T2_star = getT2Values(metabolites);

    % Build the phantom for MRSI experiements
    phantom = MRSI_build_phantom([phantom.fovY, phantom.fovX], metabolites, ...
        concentrationMap, B0, T2_star, 'use_disc', MemoryOptions.use_disc);
end



function assertIntensities(whiteMatterIntensity, greyMatterIntensity, protonDensity, skullIntensity, slice)
    % conformaty checks. Make sure pd, wm, gm are all the same size.
    if(~( ...
            isequal(size(whiteMatterIntensity), size(greyMatterIntensity)) && ...
            isequal(size(whiteMatterIntensity), size(protonDensity)) && ...
            isequal(size(whiteMatterIntensity), size(skullIntensity))) ...
            )
        error('Dimensions of MNI are not the same!!')
    end
    if(slice > min([size(whiteMatterIntensity, 3), size(greyMatterIntensity, 3), ...
            size(protonDensity, 3), size(skullIntensity, 3)]))
        error('slice index is to big')
    end
end



function phantom = calculateSpatialParameters(mniImageSize, whiteMatterCoordinates)
    %Calculate voxel size in y dimension
    firstVoxel = sub2ind(mniImageSize, 1, 1, 1);
    nextY = sub2ind(mniImageSize, 2, 1, 1);
    %Y axis is reversed (increasing y axis is moving down in the y direction)
    phantom.vox_y = whiteMatterCoordinates(1, firstVoxel) - whiteMatterCoordinates(1, nextY);
    
    %calculate voxel size in x direction
    nextX = sub2ind(mniImageSize, 1, 2, 1);
    phantom.vox_x= whiteMatterCoordinates(2, nextX) - whiteMatterCoordinates(2, firstVoxel);
    
    %calculate voxel size in z direction
    nextZ = sub2ind(mniImageSize, 1, 1, 2);
    phantom.vox_z =whiteMatterCoordinates(3, nextZ) - whiteMatterCoordinates(3, firstVoxel);

    %Calculating fov in each dimension
    phantom.fovX = phantom.vox_x * mniImageSize(2);
    phantom.fovY = phantom.vox_y * mniImageSize(1);
    phantom.fovZ = phantom.vox_z * mniImageSize(3);
end



function [whiteMatterConc, greyMatterConc] = getConcentrations
    metaboliteNames = {'PCh', 'GPC', 'Cr', 'Gln', 'Glu', 'Ins', 'NAAG', 'NAA', 'lipids', 'water'};
    %Concentrations taken from Petra J. W. Pouwels and Jens Frahm 1998. Average
    %of the 3 white matter areas.
    PchWhiteMatter = 0.8500;
    GPCWhiteMatter = 0.8500;
    CrWhiteMatter = 5.633333333333333;
    GlnWhiteMatter = 1.833333333333333;
    GluWhiteMatter = 6.566666666666666;
    InsWhiteMatter = 3.666666666666667;
    NAAGWhiteMatter = 2.266666666666667;
    NAAWhiteMatter = 7.966666666666668;
    concWhiteMatter = [PchWhiteMatter, GPCWhiteMatter, CrWhiteMatter, GlnWhiteMatter, ...
        GluWhiteMatter, InsWhiteMatter, NAAGWhiteMatter, NAAWhiteMatter, 1, 1];
    
    whiteMatterConc = containers.Map(metaboliteNames, concWhiteMatter);
    %create wm dictionary
    
    %Concentrations taken from Petra J. W. Pouwels and Jens Frahm 1998. Average
    %of the grey matter areas
    PchGreyMatter = 0.5825;
    GPCGreyMatter = 0.5825;
    CrGreyMatter = 6.700000000000000;
    GlnGreyMatter = 4.250000000000000;
    GluGreyMatter = 8.524999999999999;
    InsGreyMatter = 4.350000000000000;
    NAAGGreyMatter = 0.850000000000000;
    NAAGreyMatter = 8.450000000000000;
    
    concGreyMatter = [PchGreyMatter, GPCGreyMatter, CrGreyMatter, GlnGreyMatter, GluGreyMatter, ...
        InsGreyMatter, NAAGGreyMatter, NAAGreyMatter, 1, 1];
    greyMatterConc = containers.Map(metaboliteNames, concGreyMatter);
end


%make sure coordinate systems are the same for wm, gm, pd, and skull.
function assertCoordinates(whiteMatterCoords, greyMatterCoords, protonDensity, skullIntensity)
    if ~(isequal(whiteMatterCoords, greyMatterCoords) && ...
         isequal(whiteMatterCoords, protonDensity) && ...
         isequal(whiteMatterCoords, skullIntensity))
        error('Cordinate systems are not the same for MNI maps. Aborting!')
    end
end

% T2 values based on 3T from Stanislav Motyka et al. 2019 simulations
function T2_star = getT2Values(metabolites)
    metaboliteNames = {'NAA', 'NAAG', 'Ins', 'Cr', 'GPC', 'PCh', 'Glu', 'Gln', 'Lipids', 'H2O'};
    T2RelaxationTimes = [0.262, 0.262, 0.262, 0.150, 0.270, 0.270, 0.260, 0.260, 0.050, 1];
    t2Dictionary = containers.Map(metaboliteNames, T2RelaxationTimes);
    T2_star = zeros(length(metabolites), 1);
    for i = 1:length(metabolites)
        metaboliteLabel = metabolites(i).name;
        metaboliteName = extractMetaboliteName(metaboliteLabel);
        % get the metabolite name
        T2_star(i) = t2Dictionary(metaboliteName);
    end
    
end

function [whiteMatterIntensity, greyMatterIntensity, protonDensity, ...
                    skullIntensity, whiteMatterCoordinates] = readAtlas(slice, resolution)
    [filePath, ~, ~] = fileparts(mfilename('fullpath')); 
    if(contains(resolution, '2mm'))
        resolutionPath = '/../MNI152_2mm/';
    elseif(contains(resolution, '1mm'))
        resolutionPath = '/../MNI152_1mm/';
    else
        error('unknown resolution');
    end
    wm = spm_vol([filePath, resolutionPath, 'wmMNI152_T1_', resolution, '_brain.nii']);
    gm = spm_vol([filePath, resolutionPath, 'gmMNI152_T1_', resolution, '_brain.nii']);
    pd = spm_vol([filePath, resolutionPath, 'pdMNI152_T1_', resolution, '_brain.nii']);
    skull = spm_vol([filePath, resolutionPath, 'MNI152_T1_', resolution, '_skull.nii']);
    
    % First dimension: sagital plane (x), second dimension: coronal plane (y)
    % third dimension: transverse plane (z)
    [whiteMatterIntensity, whiteMatterCoordinates] = spm_read_vols(wm);
    [greyMatterIntensity, greyMatterCoordinates] = spm_read_vols(gm);
    [protonDensity, protonDensityCoordinates] = spm_read_vols(pd);
    [skullIntensity, skullInensityCoordinates] = spm_read_vols(skull);
    
    assertIntensities(whiteMatterIntensity, greyMatterIntensity, protonDensity, skullIntensity, slice);
    assertCoordinates(whiteMatterCoordinates, greyMatterCoordinates, protonDensityCoordinates, skullInensityCoordinates);
    % use wm intensity here beacuse we asserted earlier that they are all the
    % smae size
    
    
    % Permute to y dimension to be first. (Follows MATLAB's image dimension
    % order). Flip first dimension to have front of head first. (ie. if
    % plotted using imagesc() the front of the brain will be at the top)
    whiteMatterIntensity = permute(whiteMatterIntensity, [2,1,3]);
    greyMatterIntensity = permute(greyMatterIntensity, [2,1,3]);
    protonDensity = permute(protonDensity, [2,1,3]);
    skullIntensity = permute(skullIntensity, [2,1,3]);

    % Flip along y dimesnion to achieve right coordinate system based on matlab
    % imagaging.
    whiteMatterIntensity = flip(whiteMatterIntensity, 1);
    greyMatterIntensity = flip(greyMatterIntensity, 1);
    protonDensity = flip(protonDensity, 1);
    skullIntensity = flip(skullIntensity, 1);
end