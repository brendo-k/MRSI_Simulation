
%trajectory: trajectory array with the first dimension being excitation number
%and the second being points along resepective scan's trajectory. (ie.
%tradition MRSI's first dimension would be each excitation (in order) while the second is
%the position during the readout of that excitation. Numbers are complex to
%represent x and y positions
%
% Input:
% traj: k_space trajectory
% gMax: gradient max [mT/m] -> T/mm
function [gradientTraj, gradientTime] = MRSI_load_ktrajectory(traj, gradientMax)

    gradientMax = gradientMax/10^6; %[T/mm]

    %dimension is TR, trajectoryPoints
    kSpaceTrajectory = traj.k_trajectory; %[mm^-1]

    %initalize graident trajectory
    gradientTraj = complex(zeros(size(kSpaceTrajectory)), 0);
    %setting the first point to be complex and looping from the end has a speed
    %benefit for some reason.
    gradientTraj(1) = 0 + 1i;
    %initalize time array
    gradientTime = zeros(size(kSpaceTrajectory));

    %Loop backwards trough excitations for computing speed
    for excite = size(kSpaceTrajectory, 1):-1:1
        %Initalize the previous k space point
        previousKSpacePosition = 0;
        for k = 1:size(kSpaceTrajectory, 2)
            %Get the difference between two k space points
            kSpaceDifference = kSpaceTrajectory(excite, k) - previousKSpacePosition;

            if(k == 1)
                %calculate how to get to k_diff with the smallest amount of time
                timeStep = getFastestTime(kSpaceDifference, gradientMax);
                gradient = calculateGradient(kSpaceDifference, timeStep);
            else
                %Calculate the gradient from the difference (Numerical derrivative)
                gradient = calculateGradient(kSpaceDifference, traj.dwellTime);
                timeStep = traj.dwellTime;
            end
            if(abs(real(gradient)) >= gradientMax || abs(imag(gradient)) > gradientMax)
                error('MRSI_load_trajectory:gMaxExceeded', '%d exceeds gMax', gradient)
            end
            
            %Set k space time and gradient
            gradientTime(excite, k) = timeStep;
            gradientTraj(excite, k) = gradient;

            %set k_prev to the current k space position
            previousKSpacePosition = kSpaceTrajectory(excite, k);
        end
    end
    [firstGrad, firstTime] = recalcFirstKSpace(kSpaceTrajectory(:, 1), gradientTime(:, 1));
    gradientTime(:, 1) = firstTime;
    gradientTraj(:, 1) = firstGrad;

end

% Used to recalculate first k-space point to ensure gradients finish all at the
% same time. 
function [gradient, time] = recalcFirstKSpace(firstKSpacePoints, gradientTime)
    maxTime = max(gradientTime);
    time = repmat(maxTime, size(firstKSpacePoints));
    gradient = calculateGradient(firstKSpacePoints, maxTime);
    
end

function gradient = calculateGradient(kDifference, timeStep)
    gamma = getGamma('overTwoPi', true);
    %Calculate the gradient from the difference (Numerical derrivative)
    if(timeStep == 0)
        gradient = 0;
    else
        gradient = kDifference/(timeStep * gamma); 
    end
end

function timeStep = getFastestTime(kDifference, gradientMax)
    gamma = getGamma('overTwoPi', true);
    timeStep = abs(kDifference)/(gradientMax * gamma);
    timeStep = timeStep * 1.0001;
    timeStep = max(real(timeStep), imag(timeStep));
end