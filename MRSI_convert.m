
function out = MRSI_convert(fid, par)
    out.fids = fid;
    out.sz = size(fid);
    out.t = par.t;
    out.dwelltime = par.dwelltime;
    out.spectralwidth = par.sw;
    out.txfrq = 0;
    out.date = date;
    out.dims.t = 1;
    out.dims.x = 2;
    out.dims.y = 3;
    out.Bo = par.B0;
    out.seq = 'MRSI Simulation';
    out.te = par.te;
    out.tr = par.tr;
    out.pointsToLeftshift = 0;
    out.fovX = par.fovX;
    out.fovY = par.fovY;
    out.fovZ = par.fovZ;
    out.deltaX = par.deltaX;
    out.deltaY = par.deltaY;
    out.deltaZ = par.deltaZ;
    out.deltaK_X = par.k_delta_x;
    out.deltaK_Y = par.k_delta_y;
    out.fovK_X = par.k_fov_x;
    out.fovK_Y = par.k_fov_y;
    out.imageOrigin = [0 0 0];
    out.k_XCoordinates = par.k_min_x:k_delta_x:par.k_max_x;
    out.k_YCoordinates = par.k_min_y:k_delta_y:par.k_max_y;

    %FILLING IN THE FLAGS
    out.flags.writtentostruct=1;
    out.flags.gotparams=1;
    out.flags.leftshifted=0;
    out.flags.filtered=0;
    out.flags.zeropadded=0;
    out.flags.freqcorrected=0;
    out.flags.phasecorrected=0;
    out.flags.averaged=0;
    out.flags.addedrcvrs=0;
    out.flags.subtracted=0;
    out.flags.writtentotext=0;
    out.flags.downsampled=0;
    out.flags.spatialFT = 0;
    out.flags.spectralFT = 0;
    out.flags.coilCombined = 0;
end
