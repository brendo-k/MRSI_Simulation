
%trajectory: trajectory array with the first dimension being excitation number
%and the second being points along resepective scan's trajectory. (ie.
%tradition MRSI's first dimension would be each excitation (in order) while the second is
%the position during the readout of that excitation. Numbers are complex to
%represent x and y positions
function [gradientTraj] = MRSI_load_ktrajectory(traj, gMax) 

    dwellTime = traj.dwellTime;
    %convert to T/m
    gMax = gMax*10^-3;

    %gyromagnetic for H
    gyromagneticRatio = 42.577478518e6; %[Hz⋅T^−1]
    k_traj = traj.k_trajectory;
    
    %get the furthest k space position
    k_furthest_start = max( abs( [real(k_traj(:,1)); imag(k_traj(:,1)) ]));

    %time to ramp to furthest k space starting point.
    %TODO: add slew rate for ramping
    max_time = k_furthest_start/(gyromagneticRatio*gMax);


    %add enough space in trajectory to accountn for ramping to first position
    %gradientTraj(excite, i).time is the time to get from i-1 to k(excite,i)
    %gradientTraj(excite, i).G is the gradient applied during "time" to get to k(excite, i)
    gradientTraj = struct('G', cell(size(k_traj)), 'time', cell(size(k_traj)));

    %Calculate gradient trajectory
    for excite = 1:size(k_traj, 1)
        
        %Initalize the previous k space point
        k_prev = 0;
        for k = 1:size(k_traj, 2)
            %Get the difference between two k space points
            k_diff = k_traj(excite, k) - k_prev;

            %Set all first k space positions to the same echo time.
            if k == 1
                %Get the gradient to get to the first k position over max
                %timie
                gradient_diff = k_diff/(max_time * gyromagneticRatio);
                
                %set variable to be added to the struct
                time = max_time;
                grad_x = real(gradient_diff);
                grad_y = imag(gradient_diff);
            else
                %Calculate the gradient from the difference (Numerical derrivative)
                gradient_diff = k_diff/(dwellTime * gyromagneticRatio);
                
                %set default calculated variables
                grad_x = real(gradient_diff);
                grad_y = imag(gradient_diff);
                time_x = dwellTime;
                time_y = dwellTime;
                
                %If x gradient component is larger than gMax. 
                %Calculate how to get to same position within the contraints of gMax
                if(grad_x >= gMax || grad_x <= -gMax)
                    
                    %Get the time to get to next k posion at gMax
                    time_x = abs(grad_x)*dwellTime/gMax;
                    %set gradient to be positive or negative gmax
                    if(grad_x < 0)
                        grad_x = -gMax;
                    else
                        grad_x = gMax;
                    end
                end
                
                %If gradient is larger than gMax calculate how to get to same position with the contraints of gMax
                if(grad_y >= gMax || grad_y <= -gMax)
                    %Get the time to get to next k posion at gMax
                    time_y = abs(grad_y)*dwellTime/gMax;
                    %set grad to be positive or negative gmax
                    if(grad_y < 0)
                        grad_y = -gMax;
                    else
                        grad_y = gMax;
                    end
                end
                
                %calculate the gradient for shorter time if it was applied
                %over the longer time period.
                if(time_x < time_y)
                    grad_x = time_x*grad_x/time_y;
                elseif(time_y < time_x)
                    grad_y = time_y*grad_y/time_x;
                end
                
                %set the time to be the largest time
                time = max(time_x, time_y);
            end
            
            %Create struct from set variables
            gradient_struct.time = time;
            gradient_struct.G = grad_x + 1i*grad_y;
            
            %add created struct to array
            gradientTraj(excite, k) = gradient_struct;
            %set k_prev to the current k space position
            k_prev = k_traj(excite, k);
        end
    end
end
