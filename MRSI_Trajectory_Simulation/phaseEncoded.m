%Simulating the trjaectory of the basic cartesian MRSI trajectory
%Input: par (optional)
%           consists of fields:
%             dwellTime: dwell time in seconds
%             Fov: desired fov size in meters
%             imageSize: desired resolution. Array with 2 entries corresponding to x and y axises respectivly
%             readOutTime: readout time in the spectral dimension
%             repetitionTime: Time between each excitation    

function obj = phaseEncoded(trajectoryParameters)
arguments
    trajectoryParameters.spectralWidth (1,1) double = 2000 %[Hz]
    trajectoryParameters.Fov (2,1) double =  [200, 200]; %FoV in mm
    trajectoryParameters.imageSize (2,1) double = [8, 8]; %number of voxels in the x and y direction
    trajectoryParameters.spectralSize (1, 1) double = 512
end
    spectralWidth = trajectoryParameters.spectralWidth;
    Fov = trajectoryParameters.Fov;
    imageSize = trajectoryParameters.imageSize;
    spectralSize = trajectoryParameters.spectralSize;

    %updating local variables from parameters
    deltaFovX = Fov(1)/imageSize(1); %[mm]
    deltaFovY = Fov(2)/imageSize(2); %[mm]
    deltaKX = 1/Fov(1); %[mm^-1]
    deltaKY = 1/Fov(2); %[mm^-1]
    FovKX = 1/deltaFovX; %[mm^-1]
    FovKY = 1/deltaFovY; %[mm^-1]
    spectrallDwellTime = 1/spectralWidth; %[s]
    
    %calculating trajectory for each shot
    kSpaceX = -FovKX/2+deltaKX/2:deltaKX:FovKX/2-deltaKX/2;
    kSpaceY = -FovKY/2+deltaKY/2:deltaKY:FovKY/2-deltaKY/2;
    if(isnan(kSpaceY))
        kSpaceY = [0];
    end
    if(isnan(kSpaceX))
        kSpaceY = [0];
    end
    [meshx, meshy] = meshgrid(kSpaceX, kSpaceY);
    meshy = meshy .* 1i;
    traj = meshx + meshy;
    
    dim = size(traj);
    excite = dim(1) * dim(2);
    
    traj = reshape(traj, [excite ,1]);
    t = 0:spectrallDwellTime:spectrallDwellTime*spectralSize-spectrallDwellTime;
    
    %add plus one because the zeroth point is counted
    traj = repmat(traj, [1, spectralSize]);
    
    obj = Trajectory('cartesian', traj, imageSize, spectralSize, Fov, spectrallDwellTime, spectralWidth, t, 1);
end
    
