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
function [out, voxel_sig] = MRSI_simulate(traj, phantom, gMax, B0, MemoryOptions)
arguments
    traj (1,1) Trajectory
    phantom (:, :) struct
    gMax (1,1) double = 30
    B0 (1,1) double = 3
    MemoryOptions.use_disc (1,1) logical = 0;
end
func = tic;

if(isempty(phantom.spins{1}))
    MemoryOptions.use_disc = 1;
end


%Calculate gradient, k space, and spatical parameters
[gradient, gradientTime] = MRSI_load_ktrajectory(traj, gMax);

%create an matrix of x and y coordinates (used for speedup)
[phan_x, phan_y] = meshgrid(single(phantom.x), single(phantom.y));

%Phantom size

%Initalize array for signal readout
S = zeros(length(phantom.met), size(traj.k_trajectory, 1),  size(traj.k_trajectory,2));

S = complex(S, 0);

readout_length = size(gradient, 2);


%voxel_sig = complex(zeros(size(gradient,1), size(gradient, 2), x_phantom_size, y_phantom_size), 0);
for m = 1:length(phantom.met)
    if(MemoryOptions.use_disc)
        spins = load_spins(phantom, m);
    else
        spins = phantom.spins{m};
    end
    spins = MRSI_excite(spins, 90, 'y', 'argument_type', 'matrix', 'F', phantom.met(m).Fy);
    fprintf("simulating metabolite %s\n", phantom.met_names{m})
    scale = 2^(2-phantom.met(m).nspins);

    trc = single(phantom.met(m).Fx + 1i*phantom.met(m).Fy);
    met = phantom.met(m);
    
    spins = permute(spins, [3,4,1,2]);
    %vectorize the position dimension. Now excite spins is a stack of
    %matricies
    spins = reshape(spins, size(spins,1), size(spins,2), []);
    %excite_spins = gpuArray(excite_spins);
    idx = any(spins, [1,2]);
    spins =spins(:,:, idx);

    met_signal = complex(zeros(size(gradient, 1,2)), 0);
    traj_name = traj.name;
    ticBytes(gcp)
    parfor excite=1:size(gradient,1)
        tic
        fprintf("simulating TR number %d\n", excite)
        %permute so the spin matrix is the first two dimensions and the
        %second 2 are the voxel position dimensions

        par_spins = spins;
        for k=1:readout_length
            
            %apply gradient to x and y directions to form a gradient matrix.
            %grad_matrix has gradient strength at the x and y position;
            grad_matrix = real(gradient(excite,k))*phan_x + imag(gradient(excite,k))*phan_y + B0;
            
            %Calculate hamiltonians
            if(~(strcmp(traj_name, 'cartesian') && k ~= 1))
                [H, H_inv, is_diag] = calculate_H(met, grad_matrix, gradientTime(excite,k), B0, idx);
            end

            if(is_diag)
                par_spins = par_spins.*reshape(H, 1, size(H,1), size(H,2));
                par_spins = par_spins.*reshape(H_inv, size(H,1), 1, size(H,2));
            else
                par_spins = pagemtimes(pagemtimes(H_inv, par_spins), H);
            end

            %save to signal
            met_signal(excite, k) = readout(par_spins, trc, scale);
        end
        met_time = toc;
        fprintf("finished excitation number %d, took %f minutes\n", excite, met_time/60)
    end 
    
    tocBytes(gcp)
    t = 0:traj.dwellTime:traj.dwellTime*(size(traj.k_trajectory, 2)-1);
    S(m,:, :) = met_signal .* exp(-t/phantom.T2(m));
end
S = sum(S, 1);

S = MRSI_regrid(S, traj);
%convert to fid-a structure
out = MRSI_convert(S, traj, B0);

total = toc(func);
fprintf('total run time is %f minutes\n', total/60)
end

%FUNCTION readout: readout the signal from the spin system from all voxels
%get the readout from all voxels for metabolite m
function phantom_sig = readout(spins, Trc, scale)
traced_spins = pagemtimes(spins, Trc);
sum_sig = sum(traced_spins, 3);
phantom_sig = scale*trace(sum_sig);
end


%FUNCTION calculate_H: create hamiltonians for each voxel based on gradient strength in
%grad matrix
function [H, H_inv, is_diag] = calculate_H(met, gradients, time, B0, idx)
%get the new HAB for each spin at each voxel.
%This is derrived from ppm = ((v - v_ref)/v_ref)*10^6.
%
% (ppm*v_ref/10^6 + v_ref) = v. Where v_ref = (B0+G_x*x + G_y*y)*gamma
%
%(((ppm*v_ref/10^6 + v_ref)- v_ref_2)/v_ref_2)*10^6 = ((v-v_ref_2)/v_ref_2)*10^6 = dI_eff; where v_ref_2 = B0*gamma
%
%ppm_eff = (ppm*(B0+G*r)/10^6 + (B0+G*r) - B0)*10^6/B0
%(gamma can be factored out)

%gyromagnetic ratio
gamma=42577000;  %[Hz/T]
gradients = gradients(:);
%compute the ppm_eff from above
ppm_eff = (gradients.*met.shifts(1)/1e6 + gradients - B0)*1e6/B0;
%remove empty voxels
ppm_eff = ppm_eff(idx);
ppm_diff = ppm_eff - met.shifts(1);

%convert to rads/s
shift_rads = ppm_diff*2*pi*gamma*B0/1e6;
shift_rads = reshape(shift_rads, 1,length(shift_rads));

diag_idx = get_diag_index(size(met.Fz,1),1);
Fz = met.Fz(diag_idx);
%diag_idx = get_diage_index(size(Fz,1), 1);
Fz = repmat(Fz, [1,size(ppm_eff, 1)]); 

%element wise multiplication of shift rads. 
add_HAB = Fz.*shift_rads;

%sum along the third dimension (stack of matricies) and add HABJonly. 
%Final dimensions are (size of Hamiltonian along first dimension, size of
%Hamiltonian along second, and length of vectorized voxel positions)
new_HAB = repmat(single(met.HAB), [1,1,sum(idx)]);
diag_idx = get_diag_index(size(new_HAB, 1), size(new_HAB,3));
new_HAB(diag_idx) = new_HAB(diag_idx) + squeeze(add_HAB);
if(isdiag(new_HAB(:,:,1)))
    m = size(new_HAB,1);
    n = size(new_HAB,3);
    diag_idx = get_diag_index(m,n);

    HAB = new_HAB(diag_idx).*(1i*time);
    %initalize array for hamiltonian
    
    %get linearlized indexes of diagonals
    %extract diagonals from stack of matrices from idx
    exp_HAB = exp(HAB);
    exp_conj = conj(exp_HAB);
    H = exp_HAB;
    H_inv = exp_conj;
    is_diag = true;
    
else    
    H = complex(zeros(size(new_HAB), 'like', single(1i)));
    H(1) = 0 + 1i;
    HAB = new_HAB.*(1i*time);
    for i = size(H, 3):-1:1
        H(:,:, i) = expm(HAB(:,:,i));
    end
    H_inv = conj(H);
    is_diag = false;
end
end
