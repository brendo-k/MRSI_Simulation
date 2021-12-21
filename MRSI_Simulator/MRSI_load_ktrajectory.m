
%trajectory: trajectory array with the first dimension being excitation number
%and the second being points along resepective scan's trajectory. (ie.
%tradition MRSI's first dimension would be each excitation (in order) while the second is
%the position during the readout of that excitation. Numbers are complex to
%represent x and y positions
%
% Input:
% traj: k_space trajectory
% gMax: gradient max [mT/m] -> T/mm
function [gradientTraj, gradientTime] = MRSI_load_ktrajectory(traj, gMax)

    dwellTime = traj.dwellTime; %[s]
    gMax = gMax/(10^6); %[T/mm]

    %gyromagnetic for H
    gamma = 42577000; %[Hz⋅T^−1]
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
            k_diff = kSpaceTrajectory(excite, k) - previousKSpacePosition;

            %Calculate the gradient from the difference (Numerical derrivative)
            gradient = calculateGradient(k_diff, dwellTime, gamma);
            time = dwellTime;
            if(gradient > gMax)
                error('MRSI_load_trajectory:gMaxExceeded', '%d exceeds gMax', gradient)
            end
            

            %Create struct from set variables
            gradientTime(excite, k) = time;
            gradientTraj(excite, k) = gradient;

            %set k_prev to the current k space position
            previousKSpacePosition = kSpaceTrajectory(excite, k);
        end
    end
end

function [gradient] = calculateGradient(k_diff, dwellTime, gamma)
    %Calculate the gradient from the difference (Numerical derrivative)
    gradient = k_diff/(dwellTime * gamma);
end