
function out = MRSI_convert(fid, traj, B0)
    out.fids = fid;
    out.sz = size(fid);
    out.t = 0:1/traj.sw:1/traj.sw*(traj.imageSize(3)-1);
    out.dwelltime = traj.dwellTime;
    out.spectralwidth = traj.sw;
    out.txfrq = 0;
    out.date = date;
    out.dims.t = 1;
    out.dims.x = 2;
    out.dims.y = 3;
    out.dims.z = 0;
    out.dims.averages = 0;
    out.dims.coils = 0;
    out.Bo = B0;
    out.seq = 'MRSI Simulation';
    %TODO: Calculate  te in laod gradient
    out.te = 0;
    out.tr = traj.TR;
    out.pointsToLeftshift = 0;
    out.fovX = traj.FoV.x;
    out.fovY = traj.FoV.y;
    out.fovZ = 1;
    out.deltaX = traj.pixel_width.x; %[mm]
    out.deltaY = traj.pixel_width.y; %[mm]
    out.deltaZ = 1;
    out.imageOrigin = [0 0 0];
    out.affine_matrix = [out.deltaX, 0         , 0         , 0;...
                         0         , -out.deltaY, 0        , 0;...
                         0         , 0         , out.deltaZ, 0;...
                         0         , 0         , 0         , 1];
    out.affine_matrix = out.affine_matrix * [1, 0, 0, -out.fovX/2;...
                                             0, 1, 0, out.fovY/2;...
                                             0, 0, 1, -out.fovZ/2;...
                                             0, 0, 0 , 1];


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
