%Simulating the trjaectory of the basic cartesian MRSI trajectory
%Input: par (optional)
%           consists of fields:
%             dwellTime: dwell time in seconds
%             Fov: desired fov size in meters
%             imageSize: desired resolution. Array with 2 entries corresponding to x and y axises respectivly
%             readOutTime: readout time in the spectral dimension
%             repetitionTime: Time between each excitation    

function obj = phaseEncoded(sw, Fov, imageSize)
arguments
    sw (1,1) double = 2500 %[Hz]
    Fov (2,1) double =  [250, 250]; %FoV in mm
    imageSize (3,1) double = [8, 8, 512]; %number of voxels in the x and y direction
end

    %calculating the same image parameters as the default parameters in Rosette.m 

    %updating local variables from parameters
    deltaFovX = Fov(1)/imageSize(1); %[mm]
    deltaFovY = Fov(2)/imageSize(2); %[mm]
    deltaKX = 1/Fov(1); %[mm^-1]
    deltaKY = 1/Fov(2); %[mm^-1]
    FovKX = 1/deltaFovX; %[mm^-1]
    FovKY = 1/deltaFovY; %[mm^-1]
    dwellTime = 1/sw; %[s]
    
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
    t = 0:dwellTime: dwellTime*imageSize(3)-dwellTime;
    
    %add plus one because the zeroth point is counted
    traj = repmat(traj, [1,imageSize(3)]);
    
    obj = Trajectory('cartesian', traj, imageSize, Fov, dwellTime, sw, t, 1);
end
    
