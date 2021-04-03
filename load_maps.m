function [phantom] = load_maps()
    %Load maps
    csf = spm_vol('mni_icbm152_nlin_sym_09a/mni_icbm152_csf_tal_nlin_sym_09a.nii');
    wm = spm_vol('mni_icbm152_nlin_sym_09a/mni_icbm152_wm_tal_nlin_sym_09a.nii');
    gm = spm_vol('mni_icbm152_nlin_sym_09a/mni_icbm152_gm_tal_nlin_sym_09a.nii');

    %first dimension: sagital plane (x)
    %second dimension: coronal plane (y)
    %last dimension: transverse plane (z)
    [phantom.csf, phantom.csf_mat] = spm_read_vols(csf);
    [phantom.wm, phantom.wm_mat] = spm_read_vols(wm);
    [phantom.gm, phantom.gm_mat] = spm_read_vols(gm);

    %Calculate voxel size by taking the difference of the coordinates of the
    %adjacent voxel. All maps are the same size, thus using csf.dim to
    %calculate voxel size. If not, need to scale to same coordinates

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
end





