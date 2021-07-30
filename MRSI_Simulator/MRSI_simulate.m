% MRSI_simulate.m
% Brenden Kadota, Sunnybrook Hospital 2021
%
% This function takes in an arbitrary MRSI trajectory object and simulates
% the resulting signal from the k-space trajectory. Uses hamiltonians and
% spin systems to simulate output.
%
%
%Inputs:
% Traj: k-space trajectory object from MRSI simulations (use MRSI
% trajectory simulation package)
% gMax: Maximum gradient strength in mT/m [mT/m]
% B0: B0 magnetic field in T [T]
% T2star: T2star weighting in seconds [s]
%
% Output:
% out: FID-A MRSI object
% voxel_sig: REDUNDANT, to be used for annimations. WORK IN PROGRESS
function [out, voxel_sig] = MRSI_simulate(traj, phantom, gMax, B0, T2star)
    arguments
        traj (1,1) Trajectory
        phantom (:, :) struct 
        gMax (1,1) double = 10
        B0 (1,1) double = 3
        T2star (1,1) double = 0.1
    end
tic

%gyromagnetic ratio
gamma=42577000;  %[Hz/T]

%Calculate gradient, k space, and spatical parameters
[gradient] = MRSI_load_ktrajectory(traj, gMax);


%create an matrix of x coordinates (used for speedup)
phan_x = [phantom.x];
phan_x = reshape(phan_x, [size(phantom, 1), size(phantom, 2)]);

%create an matrix of y coordinates (used for speedup)
phan_y = [phantom.y];
phan_y = reshape(phan_y, [size(phantom, 1), size(phantom, 2)]);

%Initalize array for signal readout
S = zeros(size(traj.k_trajectory, 1), size(traj.k_trajectory,2));

S = complex(S, 0);

readout_length = size(gradient, 2);
x_phantom_size = size(phantom, 1);
y_phantom_size = size(phantom, 2);

%Now start readout:
TE = 0.001;
phantom = MRSI_excite(phantom, 90, 'y');
phantom = MRSI_evolve(phantom, TE);
phantom = MRSI_excite(phantom, 180, 'x');

%voxel_sig = complex(zeros(size(gradient,1), size(gradient, 2), x_phantom_size, y_phantom_size), 0);

parfor excite=1:size(gradient,1)
    fprintf("simulating excitation number %d\n", excite)
    
    %Excite phantom
    new_phantom = phantom;
    
    for k=1:readout_length
        %Get current gradient and k space positions
        cur_grad = gradient(excite,k).G;
        cur_time = gradient(excite,k).time;
        
        %apply gradient to x and y directions to form a gradient matrix.
        %grad_matrix has gradient strength at the x and y position;
        grad_matrix = real(cur_grad)*phan_x + imag(cur_grad)*phan_y + B0;
        
        phantom_sig = zeros(size(new_phantom));
        
        %loop through phantom
        for x = 1:x_phantom_size
            for y = 1:y_phantom_size
                
                %get the voxel's metabolites
                voxel = new_phantom(x,y).met;
                
                %loop through metabolites
                for m = 1:length(voxel)
                    
                    %get the new ppm effective of each spin. 
                    %This is derrived from ppm = ((v - v_ref)/v_ref)*10^6.
                    %
                    % (ppm*v_ref/10^6 - v_ref) = v. Where v_ref = (B0+G_x*x + G_y*y)*gamma
                    %
                    %(((ppm*v_ref/10^6 - v_ref)- v_ref_2)/v_ref_2)*10^6 = ((v-v_ref_2)/v_ref_2)*10^6 = dI_eff; where v_ref_2 = B0*gamma
                    %
                    %ppm_eff = (ppm*(B0+G*r)/10^6 + (B0+G*r) - B0)*10^6/B0
                    %(gamma can be factored out)
                    
                    dI_eff = (new_phantom(x,y).met(m).shifts*grad_matrix(x,y)/1e6 + grad_matrix(x,y) - B0)*1e6/B0;
                    
                    %get the shift in radians
                    shift_rads = dI_eff*2*pi*gamma*B0/1e6;
                    %Iz spin state
                    Iz = new_phantom(x,y).met(m).Iz;
                    %reshape shift_rads 
                    shift_rads = reshape(shift_rads, [1,1,length(shift_rads)]);
                    
                    %Calculate the HAB operator from the new frequencies
                    HAB =  squeeze(sum(Iz.*shift_rads,3));
                    %Add the j coupling to HAB operator
                    HAB = HAB + new_phantom(x,y).met(m).HABJonly;
                    
                    %create the hamiltonian
                    matrix_exp = expm(-1i*HAB*cur_time);
                    %inverse hamiltonian is just the complex conjugate
                    inverse_exp = conj(matrix_exp);
                    
                    %apply sandwich operation
                    new_phantom(x,y).d{m} = matrix_exp*new_phantom(x,y).d{m}...
                        *inverse_exp;
                   
                    %scaling factor
                    val=2^(2-new_phantom(x,y).met(m).nspins);
            
                    if k~=1
                        %save to signal
                        phantom_sig(x,y) = phantom_sig(x,y) + val*trace(new_phantom(x,y).met(m).Trc * new_phantom(x,y).d{m});
                        %voxel_sig(excite, k, x, y) = phantom_sig(x,y);
                    end
                end
            end
        end
        
        if k == 1
            %refocusing pulse for spin echo
            new_phantom = MRSI_evolve(new_phantom, TE-cur_time);
            for x = 1:x_phantom_size
                for y = 1:y_phantom_size
                    for m = 1:length(new_phantom(x,y).met)
                        %get signal after refocusing pulse
                        phantom_sig(x,y) = phantom_sig(x,y) + trace(new_phantom(x,y).met(m).Trc * new_phantom(x,y).d{m});
                        %voxel_sig(excite, k, x, y) = phantom_sig(x,y);
                    end
                end
            end
        end
        
        S(excite, k) = sum(phantom_sig, 'all');
    end
    
end

t = 0:traj.dwellTime:traj.dwellTime*(size(traj.k_trajectory, 2)-1);
%apply T2 weighting
for excite = 1:size(S, 1)
    S(excite,:) = S(excite,:).*exp(-t/T2star); 
end

S = MRSI_regrid(S, traj);
%convert to fid-a structure
out = MRSI_convert(S, traj, B0);

toc
end

