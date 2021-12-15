
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
    gamma = 42.577478518e6; %[Hz⋅T^−1]
    %dimension is TR, trajectoryPoints
    kSpaceTrajectory = traj.k_trajectory; %[mm^-1]

    firstKSpacePoint = [real(kSpaceTrajectory(:, 1)); imag(kSpaceTrajectory(:, 1))];
    %get the furthest k space position
    k_furthest_start = max(abs(firstKSpacePoint)); %[mm^-1]
    %time to get to furthest k space starting point.
    max_time = k_furthest_start/(gamma*gMax);

    %add enough space in trajectory to accountn for ramping to first position
    %gradientTraj(excite, i).time is the time to get from i-1 to k(excite,i)
    %gradientTraj(excite, i).G is the gradient applied during "time" to get to k(excite, i)
    gradientTraj = complex(zeros(size(kSpaceTrajectory)), 0);
    gradientTraj(1) = 0 + 1i;
    gradientTime = zeros(size(kSpaceTrajectory));

    %Calculate gradient trajectory
    for excite = size(kSpaceTrajectory, 1):-1:1
        %Initalize the previous k space point
        previousKSpacePosition = 0;
        for k = 1:size(kSpaceTrajectory, 2)
            %Get the difference between two k space points
            k_diff = kSpaceTrajectory(excite, k) - previousKSpacePosition;
            if k == 1
                gradient = calculateGradient(k_diff, max_time, gamma);
                time = max_time;
            else
                %Calculate the gradient from the difference (Numerical derrivative)
                gradient = calculateGradient(k_diff, dwellTime, gamma);
                time = dwellTime;
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