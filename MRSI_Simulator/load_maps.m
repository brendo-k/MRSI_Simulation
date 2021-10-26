function [phantom] = load_maps(B0, slice, mets, T2_met, T2_lip, MemoryOptions)
arguments
    B0 (1,1) double
    slice (1,1) double {mustBePositive}
    mets (:,1) struct = [];
    T2_met (1,1) double = 0.1;
    T2_lip (1,1) double = 0.001;
    MemoryOptions.use_disc (1,1) logical = 0
end
    spins = load('spinSystems.mat');
    lipids = load('Lip.mat');
    %Load maps
    proj = currentProject;
    wm = spm_vol(char(append(proj.RootFolder, '/MNI152_2mm/wmMNI152_T1_2mm_brain.nii')));
    gm = spm_vol(char(append(proj.RootFolder, '/MNI152_2mm/gmMNI152_T1_2mm_brain.nii')));
    pd = spm_vol(char(append(proj.RootFolder, '/MNI152_2mm/pdMNI152_T1_2mm_brain.nii')));
    skull = spm_vol(char(append(proj.RootFolder, '/MNI152_2mm/MNI152_T1_2mm_skull.nii')));

    %first dimension: sagital plane (x)
    %second dimension: coronal plane (y)
    %last dimension: transverse plane (z)
    [phantom.wm, phantom.wm_mat] = spm_read_vols(wm);
    [phantom.gm, phantom.gm_mat] = spm_read_vols(gm);
    [phantom.pd, phantom.pd_mat] = spm_read_vols(pd);
    [phantom.skull, phantom.skull_mat] = spm_read_vols(skull);
    
    %permute to y dimension to be first. (Follows MATLAB's image dimension
    %order). Flip first dimension to have front of head first. (ie. if
    %plotted using imagesc() the front of the brain will be at the top)
    phantom.gm = flip(permute(phantom.gm, [2,1,3]), 1);
    phantom.wm = flip(permute(phantom.wm, [2,1,3]), 1);
    phantom.pd = flip(permute(phantom.pd, [2,1,3]), 1);
    phantom.skull = flip(permute(phantom.skull, [2,1,3]), 1);
    
    if(slice > size(phantom.gm,3))
        error('slice index is to big')
    end

    %Calculate voxel size in x dimension
    first_x = sub2ind(wm.dim, 1, 1, 1);
    last_x = sub2ind(wm.dim, 2, 1, 1);
    phantom.vox_x = phantom.wm_mat(1,first_x)-phantom.wm_mat(1,last_x);

    %calculate voxel size in y direction
    first_y = sub2ind(wm.dim, 1, 1, 1);
    last_y = sub2ind(wm.dim, 1, 2, 1);
    phantom.vox_y = phantom.wm_mat(2,last_y)-phantom.wm_mat(2,first_y);

    %calculate voxel size in z direction
    first_z = sub2ind(wm.dim, 1, 1, 1);
    last_z = sub2ind(wm.dim, 1, 1, 2);
    phantom.vox_z = phantom.wm_mat(3,last_z)-phantom.wm_mat(3,first_z);

    %Calculating fov in each dimension
    phantom.fovX = phantom.vox_x*wm.dim(2);
    phantom.fovY = phantom.vox_y*wm.dim(1);
    phantom.fovZ = phantom.vox_z*wm.dim(3);
    
    %if input metabolites are empty, add all metabolites
    if(isempty(mets))
        metabolites = [spins.sysCr, spins.sysNAA, spins.sysNAAG,...
            spins.sysPCh, spins.sysGPC, spins.sysIns, spins.sysGlu,...
            spins.sysGln, spins.sysH2O, lipids.sysLip];
    else
        metabolites = mets;
    end
    
    %Concentrations taken from Petra J. W. Pouwels and Jens Frahm 1998
    %White matter concentrations
    Cr = mean([5.7 5.7 5.5]);
    NAA = mean([8.1 8.0 7.8]);
    NAAG = mean([1.5 2.7 2.6]);
    Cho = mean([1.78 1.68 1.64]);
    Ins = mean([3.8  3.1 4.1]);
    Glu = mean([7.0  6.7 6.0]);
    Gln = mean([1.8  1.5 2.2]);
    wm_conc = [Cr, NAA, NAAG, Cho, Ins, Glu, Gln, 1, 1];
    %create wm dictionary
    
    %Concentrations taken from Petra J. W. Pouwels and Jens Frahm 1998
    %grey matter concentrations
    Cr = mean([6.4 6.5 6.9 7.0]);
    NAA = mean([7.7 8.2 9.2 8.7]);
    NAAG = mean([0.7 0.5 1.4 0.8]);
    Cho = mean([1.38 1.10 0.88 1.30]);
    Ins = mean([4.3 4.3 4.1 4.7]);
    Glu = mean([8.5 8.2 8.6 8.8]);
    Gln = mean([4.4 3.8 3.9 4.9]);
    gm_conc = [Cr, NAA, NAAG, Cho, Ins, Glu, Gln, 1, 1];
    
    conc_labels = {'Cr', 'NAA', 'NAAG', 'Cho', 'Ins', 'Glu', 'Gln', 'Lipids', 'H2O'};
    
    %no we have to scale the metabolites down if there are multiple fid-a
    %structures for one metabolite. We want the sum of all spins to equal
    %the concentration
    met_names = {metabolites.name};
    %Fid-a spins have the metabolite name first then a '_' then part of the
    %molecule (ie. CH2). Just extract metabolite name into occurrences.
    occurrences = regexp(met_names, '^([a-zA-Z1-9]*)_?', 'tokens');
    %convert cell of cells to cell array of strings;
    occurrences = [occurrences{:}];
    occurrences = [occurrences{:}];
    %get number of occurences down scale concentration
    for i = 1:length(gm_conc)
        label = conc_labels{i};
        if(strcmp(label, 'Cho'))
            n_label = sum(strcmp(occurrences, 'PCh'));
            n_label = n_label + sum(strcmp(occurrences, 'GPC'));
        else
            n_label = sum(strcmp(occurrences, label));
        end
        
        %scale down the concentrations
        wm_conc(i) = wm_conc(i)/n_label;
        gm_conc(i) = gm_conc(i)/n_label;
    end

    %create met cell array which is needed for build phantom
    met_arr = cell(size(phantom.wm, [1,2]));
    %initalize to be empty
    for y = 1:size(phantom.wm,1)
        for x = 1:size(phantom.wm,2)
            counter = 1;
            for m = 1:length(metabolites)
                vox_met = metabolites(m);
                [~, name] = regexp(metabolites(m).name, '^([a-zA-Z1-9]*)_?', 'match', 'tokens');
                name = name{1}{1};
                if(strcmp(name, 'PCh') || strcmp(name, 'GPC'))
                    name = 'Cho';
                end
                idx = strcmp(conc_labels, name);
                %scale each concentration by the intensity of the mni atlas
                if(strcmp(name, 'H2O'))
                    conc = phantom.pd(y,x,slice);
                elseif(strcmp(name, 'Lipids'))
                    conc = phantom.skull(y,x,slice);
                else
                    conc = wm_conc(idx)*phantom.wm(y,x,slice) + gm_conc(idx)*phantom.gm(y,x,slice);
                end
                vox_met.conc = conc;

                %add to met array
                met_arr{y,x}(counter) = vox_met;
                counter = counter + 1;
                
            end
        end
    end
    


    T2_star = zeros(length(metabolites), 1);
    for m = 1:length(metabolites)
        [~, name] = regexp(metabolites(m).name, '^([a-zA-Z1-9]*)_?', 'match', 'tokens');
        name = name{1}{1};
        if(strcmp(name, 'Lipids'))
            T2_star(m) = T2_lip;
        else
            T2_star(m) = T2_met;
        end
    end

    
    phantom = MRSI_build_phantom([phantom.fovY, phantom.fovX], met_arr, B0, T2_star,...
        'use_disc', MemoryOptions.use_disc);
end
