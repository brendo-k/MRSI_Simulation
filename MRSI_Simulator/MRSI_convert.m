
function out = MRSI_convert(fid, traj, B0)
    out.fids = fid;
    out.sz = size(fid);
    out.t = 0:1/traj.sw:1/traj.sw*(traj.imageSize(3)-1);
    out.dwelltime = traj.dwellTime;
    out.spectralwidth = traj.sw;
    out.txfrq = 0;
    out.date = date;
    out.dims.t = 3;
    out.dims.x = 1;
    out.dims.y = 2;
    out.dims.coils = 0;
    out.Bo = B0;
    out.seq = 'MRSI Simulation';
    %TODO: Calculate  te in laod gradient
    out.te = 0;
    out.tr = traj.repetitionTime;
    out.pointsToLeftshift = 0;
    out.fovX = traj.FoV.x*1000;
    out.fovY = traj.FoV.y*1000;
    out.fovZ = 1;
    out.deltaX = traj.pixel_width.x*1000; %converting [m] to [mm]
    out.deltaY = traj.pixel_width.y*1000; %converting [m] to [mm]
    out.deltaZ = 1;
    out.deltaK_X = traj.delta_K.x/1000; %converting [m^-1] to [mm^-1]
    out.deltaK_Y = traj.delta_K.y/1000; %converting [m^-1] to [mm^-1]
    out.fovK_X = traj.FovK.x/1000;
    out.fovK_Y = traj.FovK.y/1000;
    out.imageOrigin = [0 0 0];
    out.k_XCoordinates = -out.fovK_X/2+out.deltaK_X/2:out.deltaK_X:out.fovK_X/2-out.deltaK_X/2;
    out.k_YCoordinates = -out.fovK_Y/2+out.deltaK_Y/2:out.deltaK_Y:out.fovK_Y/2-out.deltaK_Y/2;

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
