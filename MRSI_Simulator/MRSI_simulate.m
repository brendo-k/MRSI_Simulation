% MRSI_simulate.m
% Brenden Kadota, Sunnybrook Hospital 2021
%
% This function takes in an arbitrary MRSI trajectory object and simulates
% the resulting signal from the k-space trajectory. Uses hamiltonians and
% spin systems to simulate output.
%
%_trajkkjkkkj:w;j;kj;kj
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
    gMax (1,1) double = 30
    B0 (1,1) double = 3
    T2star (1,1) double = 0.1
end
tic

%gyromagnetic ratio
gamma=42577000;  %[Hz/T]

%Calculate gradient, k space, and spatical parameters
[gradient] = MRSI_load_ktrajectory(traj, gMax);

%create an matrix of x and y coordinates (used for speedup)
[phan_x, phan_y] = meshgrid(phantom.x, phantom.y);

%Phantom size
phan_size = size(phan_x);

%Initalize array for signal readout
S = zeros(size(traj.k_trajectory, 1), size(traj.k_trajectory,2));

S = complex(S, 0);

readout_length = size(gradient, 2);

%Now start readout:
TE = 0.001;
phantom = MRSI_excite(phantom, 90, 'y');
phantom = MRSI_evolve(phantom, TE);
phantom = MRSI_excite(phantom, 180, 'x');

%voxel_sig = complex(zeros(size(gradient,1), size(gradient, 2), x_phantom_size, y_phantom_size), 0);

for excite=1:size(gradient,1)
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
        phantom_sig = zeros(phan_size);
        
        %loop through phantom
        for m = 1:length(phantom.spins)
            if k~=1

                %save to signal
                phantom_sig(y,x) = phantom_sig(y,x) + readout(m, new_phantom);

            end
            
            %get the metabolites
            met = new_phantom.met(m);
            
            %get the new ppm effective of each spin.
            %This is derrived from ppm = ((v - v_ref)/v_ref)*10^6.
            %
            % (ppm*v_ref/10^6 + v_ref) = v. Where v_ref = (B0+G_x*x + G_y*y)*gamma
            %
            %(((ppm*v_ref/10^6 + v_ref)- v_ref_2)/v_ref_2)*10^6 = ((v-v_ref_2)/v_ref_2)*10^6 = dI_eff; where v_ref_2 = B0*gamma
            %
            %ppm_eff = (ppm*(B0+G*r)/10^6 + (B0+G*r) - B0)*10^6/B0
            %(gamma can be factored out)
            gradients = repmat(grad_matrix, [1,1,length(met.shifts)]);
            shifts = reshape(met.shifts, [1,1,length(met.shifts)]);
            total_shifts = (gradients.*shifts/1e6 + gradients - 3)*1e6/B0;
            
            %vectorize dI_eff
            total_shifts = reshape(total_shifts, phan_size(1)*phan_size(2), length(met.shifts));
            
            %get the shift in radians
            shift_rads = total_shifts*2*pi*gamma*B0/1e6;
            %Iz spin state
            Iz = met.Iz;
            HAB = zeros(phan_size(1)*phan_size(2),size(Iz, 1), size(Iz, 2));
            shift_rads = reshape(shift_rads, phan_size(1)*phan_size(2), 1, []);
            for i = size(shift_rads, 1)
                HAB(i,:,:) = squeeze(sum(Iz.*shift_rads(i,1,:),3)) + met.HABJonly;
            end
            
            for i = size(shift_rads,1)
                matrix_exp = expm(-1i*squeeze(HAB(i,:,:))*cur_time);
                %inverse hamiltonian is just the complex conjugate
                inverse_exp = conj(matrix_exp);
                [y,x] = ind2sub(phan_size, i);
                new_phantom.spins{m}(y,x,:,:) = matrix_exp*squeeze(new_phantom.spins{m}(y,x,:,:))...
                    *inverse_exp;
                
            end
        end
        if(k == 1)
            new_phantom = MRSI_evolve(new_phantom, TE-cur_time);
            for m = 1:length(phantom.spins)
                phantom_sig(y,x) = phantom_sig(y,x) + readout(m,new_phantom);
            end
        end
    end
    
    
    S(excite, k) = sum(phantom_sig, 'all');
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

%get the readout from all voxels for metabolite m
function phantom_sig = readout(m, phantom)
for y = 1:length(phantom.y)
    for x = 1:length(phantom.x)
        scale=2^(2-phantom.met(m).nspins);
        phantom_sig = scale*trace(new_phantom.met(m).Trc * squeeze(new_phantom.spins{m}(y,x,:,:)));
    end
end
end
