function [phantom] = load_maps(B0, slice)
    spins = load('spinSystems.mat');
    %Load maps
    csf = spm_vol('../mni_icbm152_nlin_sym_09a/mni_icbm152_csf_tal_nlin_sym_09a.nii');
    wm = spm_vol('../mni_icbm152_nlin_sym_09a/mni_icbm152_wm_tal_nlin_sym_09a.nii');
    gm = spm_vol('../mni_icbm152_nlin_sym_09a/mni_icbm152_gm_tal_nlin_sym_09a.nii');

    %first dimension: sagital plane (x)
    %second dimension: coronal plane (y)
    %last dimension: transverse plane (z)
    [phantom.csf, phantom.csf_mat] = spm_read_vols(csf);
    [phantom.wm, phantom.wm_mat] = spm_read_vols(wm);
    [phantom.gm, phantom.gm_mat] = spm_read_vols(gm);
    %Calculate voxel size by taking the difference of the coordinates of the
    %adjacent voxel. All maps are the same size, thus using csf.dim to
    %calculate voxel size. If not, need to scale to same coordinates
    phantom.csf = flip(permute(phantom.csf, [2,1,3]), 1);
    phantom.gm = flip(permute(phantom.gm, [2,1,3]), 1);
    phantom.wm = flip(permute(phantom.wm, [2,1,3]), 1);

    %Calculate voxel size in x dimension
    first_x = sub2ind(csf.dim, 1, 1, 1);
    last_x = sub2ind(csf.dim, 2, 1, 1);
    phantom.vox_x = phantom.csf_mat(1,last_x)-phantom.csf_mat(1,first_x);

    %calculate voxel size in y direction
    first_y = sub2ind(csf.dim, 1, 1, 1);
    last_y = sub2ind(csf.dim, 1, 2, 1);
    phantom.vox_y = phantom.csf_mat(2,last_y)-phantom.csf_mat(2,first_y);

    %calculate voxel size in z direction
    first_z = sub2ind(csf.dim, 1, 1, 1);
    last_z = sub2ind(csf.dim, 1, 1, 2);
    phantom.vox_z = phantom.csf_mat(3,last_z)-phantom.csf_mat(3,first_z);

    %Calculating fov in each dimension
    phantom.fovX = phantom.vox_x*csf.dim(1);
    phantom.fovY = phantom.vox_y*csf.dim(2);
    phantom.fovZ = phantom.vox_z*csf.dim(3);
    

    
    metabolites = [spins.sysCr, spins.sysNAA, spins.sysNAAG, spins.sysPCh,...
                    spins.sysGPC, spins.sysIns, spins.sysGlu, spins.sysGln, spins.sysH2O];
    wm_met = metabolites;
    gm_met = metabolites;
    
    %Concentrations taken from Petra J. W. Pouwels and Jens Frahm 1998
    Cr = mean([5.7 5.7 5.5]);
    NAA = mean([8.1 8.0 7.8]);
    NAAG = mean([1.5 2.7 2.6]);
    Cho = mean([1.78 1.68 1.64]);
    Ins = mean([3.8  3.1 4.1]);
    Glu = mean([7.0  6.7 6.0]);
    Gln = mean([1.8  1.5 2.2]);
    wm_conc = [Cr, NAA, NAAG, Cho, Ins, Glu, Gln];
    max_conc_wm = max(wm_conc);
    wm_conc = wm_conc/max_conc_wm;
    wm_conc(end + 1) = 1;
    for i = 1:length(wm_met)
        if(i ==1)
            c = 1;
        elseif(i >=2 && i <= 3)
            c = 2;
        elseif(i >= 4 && i <= 6)
            c = 3;
        elseif(i >= 7 && i <= 11)
            c = 4;
        elseif(i == 12)
            c = 5;
        elseif(i == 13)
            c = 6;
        elseif(i == 14 || i == 15)
            c = 7;
        else
            c = 8;
        end
        wm_met(i).scaleFactor = wm_met(i).scaleFactor * wm_conc(c);  
    end
    
    Cr = mean([6.4 6.5 6.9 7.0]);
    NAA = mean([7.7 8.2 9.2 8.7]);
    NAAG = mean([0.7 0.5 1.4 0.8]);
    Cho = mean([1.38 1.10 0.88 1.30]);
    Ins = mean([4.3 4.3 4.1 4.7]);
    Glu = mean([8.5 8.2 8.6 8.8]);
    Gln = mean([4.4 3.8 3.9 4.9]);

    gm_conc = [Cr, NAA, NAAG, Cho, Ins, Glu, Gln];
    max_conc_gm = max(gm_conc);
    gm_conc = gm_conc/max_conc_gm;
    gm_conc(end + 1) = 1;
    for i = 1:length(gm_met)
        if(i ==1)
            c = 1;
        elseif(i >=2 && i <= 3)
            c = 2;
        elseif(i >= 4 && i <= 6)
            c = 3;
        elseif(i >= 7 && i <= 11)
            c = 4;
        elseif(i == 12)
            c = 5;
        elseif(i == 13)
            c = 6;
        elseif(i == 14 || i == 15)
            c = 7;
        else
            c = 8;
        end
        gm_met(i).scaleFactor = gm_met(i).scaleFactor * gm_conc(c);  
    end

    phantom.wm = phantom.wm > 0.5;
    phantom.gm = phantom.gm > 0.5;

    met_arr = cell(size(phantom.wm, 1), size(phantom.wm, 2));
    met_arr(phantom.wm(:,:,slice) == 1) = {wm_met};
    met_arr(phantom.gm(:,:,slice) == 1) = {gm_met};
    phantom = MRSI_build_phantom_mem([phantom.fovY, phantom.fovX], met_arr, B0);
end
