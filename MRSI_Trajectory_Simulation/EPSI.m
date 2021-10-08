%Simulating the trjaectory of an EPSI trajectory. Only simulating EPSI
%trajectory that aquires in one direction (no aquasition during rewind
%gradients).
%Input: par (optional)
%           consists of fields:
%               dwellTime: dwell time in seconds
%               Fov: desired fov size in meters
%               imageSize: desired resolution. Array with 3 entries corresponding to x, y,and spectral axises respectivly
%               readOutTime: readout time in the spectral dimension
%               repetitionTime: Time between each excitation
%
%Output:    Traj: 2d matrix of the resulting K-space trajectory from the input
%              params. Dimensions are excitation number and readout for the frist and
%              second dimensions respectivly.
%           scanTime: Resulting scan time from trajectory. In seconds
%           par: paramaters (same fields as above) used for EPSI simulation   

function [obj] = EPSI(par)

    if(nargin < 1)
        par.sw = 2500; %[Hz]
        par.Fov = [200, 200]; %FoV in m in the x and y directions
        par.imageSize = [16, 16, 1024]; %voxels resolution in the x and y and spectral dimension
    end
    gamma = 42.577478518e6;
    %updating local variables from parameters
    Fov = par.Fov;
    deltaFovX = Fov(1)/par.imageSize(1); %[mm]
    deltaFovY = Fov(2)/par.imageSize(2); %[mm]
    deltaKX = 1/Fov(1); %[1/mm]
    deltaKY = 1/Fov(2); %[1/mm]
    FovKX = 1/deltaFovX; %[1/mm]
    FovKY = 1/deltaFovY; %[1/mm]
    dwellTime = (par.imageSize*2)/par.sw; %[s]
    
    %calculating trajectory for each shot
    kSpaceX = -FovKX/2 + deltaKX/2:deltaKX:FovKX/2 - deltaKX/2;
    kSpaceY = -FovKY/2 + deltaKY/2:deltaKY:FovKY/2 - deltaKY/2;
    
    %initalize empty array for trajectory 
    traj = zeros(par.imageSize(2), par.imageSize(1)*par.imageSize(3));
    
    %readout along the first dimension
    for j = 1:par.imageSize(2)
        %multiple reads of the x axis
        traj(j,:) = repmat(kSpaceX, [1,par.imageSize(3)]);
        
        %add the corresponding y position
        traj(j,:) = traj(j,:) + 1i*kSpaceY(j);
    end
    
    %Calculate scan time.
    %(dwelltime * imageSize * 2) is the time to get back to the same position
    %in k space which needs to be repeated equivalently to spectral points.
    obj = Trajectory(traj, par.imageSize, par.Fov, dwellTime, ["x", "y"], sw);
end
    
