%MRSI_regrid.m
%Regrids simulted signal to a cartesian axis if aquired using epsi or
%phaseEncoding.

function signalRemapped = MRSI_regrid(originalSignal, traj)
arguments
    originalSignal (:,:, :) double
    traj (1,1) Trajectory
end


if(contains(traj.name, 'cart', 'IgnoreCase', 1))
    %cartesian MRSI
    signalRemapped = zeros([traj.imageSize, traj.spectralLength]);
    signalRemapped = permute(signalRemapped, [3,1,2]);
    
    for readout = 1:size(originalSignal, 2)
        [y, x] = ind2sub(traj.imageSize(1:2), readout);
        signalRemapped(:, y, x) = originalSignal(1, readout, :);
    end
elseif(contains(traj.name, 'epsi', 'IgnoreCase', true))
    signalRemapped = zeros([traj.imageSize, traj.spectralLength]);
    signalRemapped = permute(signalRemapped, [3,1,2]);

    for iPoint = 1:traj.spectralLength
        signalStart = (traj.spatialPoints * (iPoint - 1)) + 1;
        signalEnd = (traj.spatialPoints * iPoint);
        signalRemapped(iPoint, :, :) = originalSignal(1, :, signalStart:signalEnd);
    end
else
    signalRemapped = permute(originalSignal, [3, 1, 2]);
end

    
end
