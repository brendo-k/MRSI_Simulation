% PSF.m
% Jamie Near, McGill University 2019.
% 
% USAGE:
% PSF(inTrajectory)
% 
% DESCRIPTION:
% This function takes in a trajectory strucutre and displays a graph of the
% point spread function. Furhter, it returns values for full width half
% maximum and signal leakage. The singal leakage is the sum of absolute
% value of the signal outside of the main signal peak.
%
% INPUTS:
% TrajectoryStruct = FID-A Trajectory structure from CreateTrajectory.m
%
% OUTPUTS:
% lineWidth   = the full width half maximum
% leakage     = leakage of point spread function signal outside of the main
%               peak



function [lineWidth,leakage] = PSF(trajectory)
    k_traj = trajectory.k_trajectory(:, 1:trajectory.spatialPoints);
    k_space_trajectory = [real(k_traj(:)), imag(k_traj(:))];
    
    x_mesh = trajectory.spacial_coordinates.x(1):0.001:trajectory.spacial_coordinates.x(end);
    y_mesh = trajectory.spacial_coordinates.y(1):0.001:trajectory.spacial_coordinates.x(end);
    [x,y] = meshgrid(x_mesh, y_mesh);
    pointSpreadCoordinates(:,:) = [x(:), y(:)];
    
    sftOperator = sft2_Operator(k_space_trajectory, pointSpreadCoordinates, 0);
    kSpace = ones(1,size(k_space_trajectory,1));
    pointSpreadFunction = sftOperator*kSpace(:);
    pointSpreadFunction = reshape(pointSpreadFunction, [size(x,1), size(x,2)]);
    peak = max(pointSpreadFunction, [], 'all');
    pointSpreadFunction = pointSpreadFunction./peak;

    figure;
    subplot(2,2,1)
    if(size(k_traj, 2) ~= 1)
        hold on
        for i = 1:size(k_traj, 1)
            plot(real(k_traj(i,:)), imag(k_traj(i,:))), title('kSpace trajectory');
        end
        hold off
    else
        scatter(real(k_traj),imag(k_traj));
    end
   
    
    subplot(2,2,2), mesh(x_mesh, y_mesh, ...
    real(pointSpreadFunction)), axis square, title('Point Spread Function')
    
    pointSpread2D = real(pointSpreadFunction(:,(floor(end/2)+1)));
    subplot(2,2,3), plot(x_mesh, pointSpread2D), title('2D point spread through the middle');
    
    lineWidth = fwhm(x_mesh, pointSpread2D);
    
    
    middle = find(pointSpread2D == max(pointSpread2D));
    startOfLeakage = find(pointSpread2D(middle:end) < 0, 1);
    leakage = sum(sqrt(pointSpread2D(startOfLeakage:end).^2));
    
    
end