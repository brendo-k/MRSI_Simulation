
%trajectory: trajectory array with the first dimension being excitation number
%and the second being points along resepective scan's trajectory. (ie.
%tradition MRSI's first dimension would be each excitation (in order) while the second is
%the position during the readout of that excitation. Numbers are complex to
%represent x and y positions
function [trajectory, grad, par] = load_trajectory(trajectory, dwellTime, gMax) 

    %convert to T/m
    gMax = gMax*10^-3;

    %gyromagnetic for H
    gyromagneticRatio = 42.577478518e6; 
    
    %get the furthest k space position
    k_furthest_start = max(abs([real(trajectory(:,1)); imag(trajectory(:,1))]));

    %time to ramp to furthest k space starting point.
    %TODO: add slew rate for ramping
    max_time = k_furthest_start/(gyromagneticRatio*gMax);


    %add enough space in trajectory to accountn for ramping to first position
    gradientTraj = struct('G', cell(size(trajectory)), 'time', cell(size(trajectory)));

    %Calculate gradient trajectory
    for excite = 1:size(trajectory, 1)
%         %First k space position
%         first_k = trajectory(excite, 1);
% 
%         %get ramping time needed for k space point
%         g_x = real(first_k)/(gyromagneticRatio*max_time);
%         g_y = imag(first_k)/(gyromagneticRatio*max_time);
%         G = g_x + 1i*g_y;
%         
%         gradient_struct.time = max_time;
%         gradient_struct.G = G;
% 
%         %append ramping points to the first points of gradient trajectory
%         gradientTraj(excite,1) = gradient_struct;   
        %Calculate remaining points with derivative
        past = 0;
        for k = 1:length(trajectory(excite,:))
            %Get the difference between two k space points
            k_diff = trajectory(excite, k) - past;
            %Calculate the gradient from the difference (Numerical derrivative)
            gradient_diff = k_diff/(dwellTime * gyromagneticRatio);

            %If gradient is larger than gMax calculate how to get to same position with the contraints of gMax
            if(abs(real(gradient_diff)) > gMax)
                time_x = abs(real(gradient_diff))*dwellTime/gMax;
                if(real(gradient_diff) < 0)
                    grad_x = -gMax;
                else
                    grad_x = gMax;
                end
            else
                time_x = dwellTime;
                grad_x = real(gradient_diff);
            end

            %If gradient is larger than gMax calculate how to get to same position with the contraints of gMax
            if(abs(imag(gradient_diff)) > gMax)
                time_y = abs(imag(gradient_diff))*dwellTime/gMax;
                if(imag(gradient_diff) < 0)
                    grad_y = -gMax;
                else
                    grad_y = gMax;
                end
            else
                time_y = dwellTime;
                grad_y = imag(gradient_diff);
            end

            if(time_x < time_y)
                grad_x = time_x*grad_x/time_y;
            elseif(time_y < time_x)
                grad_y = time_y*grad_y/time_x;
            end
            
            gradient_struct.time = max(time_x, time_y);
            gradient_struct.G = grad_x + 1i*grad_y;

            gradientTraj(excite, k) = gradient_struct;
            past = trajectory(excite, k);
        end
    end
    grad = gradientTraj;

    %coordinates of the spacial points
    spacial_points = [];

    %loop through each excitation
    for excite = 1:size(trajectory, 1)
        %first point of the excitation
        first = trajectory(excite,1);

        %add first point to spacial points
        spacial_points(end+1) = first;

        %when trajectory(i,excite) == first we have made one cycle in k space.
        %stop looping when first is reached as remaning points are the same.
        for i=2:length(trajectory(excite,:))
            if (trajectory(excite, i) ~= first)
                spacial_points(end+1) = trajectory(excite, i);
            else
                break;
            end
        end
    end

    %sort by the y axis  
    [~,index] = sort(imag(spacial_points));
    sorted_spacial = spacial_points(index);

    %get the first point's y coordinate
    y_coord = imag(sorted_spacial(1));

    %loop through the spacial points 
    for i = 1:length(sorted_spacial)
        %when y coord is no longer the same we have the length of the x axis
        if y_coord ~= imag(sorted_spacial(i))
            row_length = i-1;
            break;
        end
    end
    %Assuming cartesian trajectories for now

    %TODO: Figure out about sampling for non cartesian
    column_length = length(sorted_spacial)/row_length;

    %set nuber of kx and ky points
    k_nx = column_length;
    k_ny = row_length;

    %Get k_max and k_min
    k_max_x = max(real(trajectory), [], 'all');
    k_max_y = max(imag(trajectory), [], 'all');
    k_min_x = min(real(trajectory), [], 'all');
    k_min_y = min(imag(trajectory), [], 'all');

    %Calculate fov (subtract k_min because negative)
    k_fov_x = k_max_x - k_min_x;
    k_fov_y = k_max_y - k_min_y;

    %calculate delta k 
    k_delta_x = k_fov_x/k_nx;
    k_delta_y = k_fov_y/k_ny;
    
    %calculate fov
    fovX = 1/k_delta_x;
    fovY = 1/k_delta_y;

    %calculate delta fov
    delta_x = 1/k_fov_x;
    delta_y = 1/k_fov_y;
    
    %update par variable
    par.k_nx = k_nx;
    par.k_ny = k_ny;
    par.k_max_x = k_max_x;
    par.k_max_y = k_max_y;
    par.k_min_x = k_min_x;
    par.k_min_y = k_min_y;
    par.k_fov_x = k_fov_x;
    par.k_fov_y = k_fov_y;
    par.k_delta_x = k_delta_x;
    par.k_delta_y = k_delta_y;
    par.fovX = fovX;
    par.fovY = fovY;
    par.delta_x = delta_x;
    par.delta_y = delta_y;
end
