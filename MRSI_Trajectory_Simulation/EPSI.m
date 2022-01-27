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

function [trajectoryObject] = EPSI(spectralWidth, Fov, imageSize, spectralSize)
arguments
    spectralWidth (1, 1) double = 2000
    Fov (2, 1) double = [200, 200]   % [y, x]
    imageSize (2, 1) double = [8, 8] % [y, x]
    spectralSize (1, 1) double = 512;
end
    numY = imageSize(1);
    numX = imageSize(2);
    fovY = Fov(1);
    fovX = Fov(2);
    
    %updating local variables from parameters
    deltaFovX = fovX/numX; %[mm]
    deltaFovY = fovY/numY; %[mm]
    deltaKX = 1/fovY; %[1/mm]
    deltaKY = 1/fovX; %[1/mm]
    FovKX = 1/deltaFovX; %[1/mm]
    FovKY = 1/deltaFovY; %[1/mm]
    
    spectralDwellTime = 1/spectralWidth;
    adcDwellTime = spectralDwellTime/((numX - 1)*2); %[s]
    
    %calculating trajectory for each shot
    kSpaceX = createCoordinates(FovKX/2, deltaKX);
    kSpaceY = createCoordinates(FovKY/2, deltaKY);
    
    %initalize empty array for trajectory 
    traj = zeros(numY, numX * spectralSize);
    
    %readout along the first dimension
    for iTrajectory = 1:numY
        %multiple reads of the x axis
        traj(iTrajectory, :) = repmat(kSpaceX, [1, spectralSize]);
        
        %add the corresponding y position
        traj(iTrajectory,:) = traj(iTrajectory,:) + 1i*kSpaceY(iTrajectory);
    end
    crossingTimeVector = 0:adcDwellTime:numX*adcDwellTime - adcDwellTime;
    timeVector = zeros(numX * spectralSize, 1);
    
    for iSpectralPoint = 1:spectralSize
        indexStart = (iSpectralPoint - 1) * numX + 1;
        indexEnd = iSpectralPoint * numX;

        nextTimeVector = (iSpectralPoint - 1) * spectralDwellTime + crossingTimeVector;
        timeVector(indexStart:indexEnd) = nextTimeVector;

    end
    %Calculate scan time.
    %(dwelltime * imageSize * 2) is the time to get back to the same position
    %in k space which needs to be repeated equivalently to spectral points.
    trajectoryObject = Trajectory('EPSI', traj, imageSize, spectralSize, Fov, adcDwellTime, spectralWidth, timeVector, numX);
end
    
