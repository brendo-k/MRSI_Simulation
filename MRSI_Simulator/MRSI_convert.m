
function out = MRSI_convert(data, traj, B0)
    out.data = data;
    out.sz = size(data);
    out.spectralTime = 0:1/traj.sw:1/traj.sw*(traj.imageSize(3)-1);
    out.spectralDwellTime = 1/traj.sw;
    out.spectralWidth = traj.sw;
    out.adcTime = traj.t;
    out.adcDwellTime = traj.dwellTime;
    out.txfrq = 0;
    out.date = date;
    out.dims.t = 1;
    out.dims.kx = 2;
    out.dims.ky = 3;
    out.dims.z = 0;
    out.dims.x = 0;
    out.dims.y = 0;
    out.dims.kz = 0;
    out.dims.extras = 0;
    out.dims.averages = 0;
    out.dims.subSpecs = 0;
    out.dims.coils = 0;
    out.Bo = B0;
    out.seq = 'MRSI Simulation';
    %TODO: Calculate  te in laod gradient
    out.te = 0;
    out.tr = traj.TR;
    out.pointsToLeftshift = 0;
    out.fov.x = traj.FoV.x;
    out.fov.y = traj.FoV.y;
    out.fov.z = 1;
    out.voxelSize.x = traj.pixel_width.x; %[mm]
    out.voxelSize.y = traj.pixel_width.y; %[mm]
    out.voxelSize.z = 1;
    out.imageOrigin = [0 0 0];

    
    out.affine_matrix = [1, 0, 0, -out.fov.x/2;...
                                             0, 1, 0, out.fov.y/2;...
                                             0, 0, 1, -out.fov.z/2;...
                                             0, 0, 0 , 1];
    out.affine_matrix = out.affine_matrix * [out.voxelSize.x, 0         , 0         , 0;...
                         0         , -out.voxelSize.y, 0        , 0;...
                         0         , 0         , out.voxelSize.z, 0;...
                         0         , 0         , 0         , 1];
    fovX = getFov(out, 'x');
    voxSizeX = getVoxSize(out, 'x');
    xCoordinates = -fovX/2 + voxSizeX/2 : voxSizeX : fovX/2 - voxSizeX/2;
    xCoordinates = xCoordinates - getImageOrigin(out, 'x');
    
    fovY = getFov(out, 'y');
    voxSizeY = getVoxSize(out, 'y');
    yCoordinates = -fovY/2 + voxSizeY/2 : voxSizeY : fovY/2 - voxSizeY/2;
    yCoordinates = yCoordinates - getImageOrigin(out, 'y');

    out.coordinates.x = xCoordinates;
    out.coordinates.y = yCoordinates;


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
    out.flags.addedrcvrs = 1;
    out.flags.isFourSteps = 0;
end
