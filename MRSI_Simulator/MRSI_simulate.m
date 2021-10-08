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
function [out, voxel_sig] = MRSI_simulate(traj, phantom, gMax, B0, T2star, MemoryOptions)
arguments
    traj (1,1) Trajectory
    phantom (:, :) struct
    gMax (1,1) double = 30
    B0 (1,1) double = 3
    T2star (1,1) double = 0.1
    MemoryOptions.use_disc (1,1) logical = 0;
end
func = tic;



%Calculate gradient, k space, and spatical parameters
[gradient] = MRSI_load_ktrajectory(traj, gMax);

%create an matrix of x and y coordinates (used for speedup)
[phan_x, phan_y] = meshgrid(single(phantom.x), single(phantom.y));

%Phantom size

%Initalize array for signal readout
S = zeros(size(traj.k_trajectory, 1), size(traj.k_trajectory,2), length(phantom.met));

S = complex(S, 0);

readout_length = size(gradient, 2);

%Now start readout:
TE = 0.001;


%voxel_sig = complex(zeros(size(gradient,1), size(gradient, 2), x_phantom_size, y_phantom_size), 0);
for m = 1:length(phantom.met)
    if(MemoryOptions.use_disc)
        spins = load_spins(phantom.file{m}, phantom, m);
    else
        spins = phantom.spins{m};
    end
    spins = MRSI_excite(spins, 90, 'y', 'argument_type', 'matrix', 'F', phantom.met(m).Fy);
    spins = MRSI_evolve(spins, TE, 'argument_type', 'matrix', 'HAB', phantom.met(m).HAB);
    spins = MRSI_excite(spins, 180, 'x', 'argument_type', 'matrix', 'F', phantom.met(m).Fx);
    fprintf("simulating metabolite %s\n", phantom.met_names{m})
    scale = 2^(2-phantom.met(m).nspins);
    trc = phantom.met(m).Fx + 1i*phantom.met(m).Fy;
    for excite=1:size(gradient,1)
        tic
        fprintf("simulating excitation number %d\n", excite)
        %permute so the spin matrix is the first two dimensions and the
        %second 2 are the voxel position dimensions
        excite_spins = permute(spins, [3,4,1,2]);
        %vectorize the position dimension. Now excite spins is a stack of
        %matricies
        excite_spins = reshape(excite_spins, size(excite_spins,1), size(excite_spins,2), []);
        %excite_spins = gpuArray(excite_spins);
        idx = any(excite_spins, [1,2]);
        for k=1:readout_length
            %if k is 1 (ie starting) don't aquire signal and wait for spin
            %echo
            if k~=1
                %save to signal
                S(excite, k, m) = readout(excite_spins, trc, scale, idx);
            end
            
            %Get current gradient and k space positions
            cur_grad = gradient(excite,k).G;
            cur_time = gradient(excite,k).time;
            
            %apply gradient to x and y directions to form a gradient matrix.
            %grad_matrix has gradient strength at the x and y position;
            grad_matrix = real(cur_grad)*phan_x + imag(cur_grad)*phan_y + B0;
            
            %Calculate hamiltonians
            [H, H_inv] = calculate_H(phantom, grad_matrix, m, cur_time, B0, idx);
            
            temp = excite_spins(:,:,idx);
            %apply hamiltonians across all voxels
            excite_spins(:,:,idx) = pagemtimes(pagemtimes(H_inv, temp), H);
            
            %aquire signal after evolution period of spin echo.
            if(k == 1)
                excite_spins = reshape(excite_spins, size(spins,3), size(spins,4), size(spins,1), size(spins,2));
                excite_spins = permute(excite_spins, [3,4,1,2]);
                
                excite_spins = MRSI_evolve(excite_spins, TE-cur_time, 'argument_type', 'matrix', 'HAB', phantom.met(m).HAB);
                excite_spins = permute(excite_spins, [3,4,1,2]);
                excite_spins = reshape(excite_spins, size(spins,3), size(spins,4), []);
                S(excite, k, m) = readout(excite_spins, trc, scale, idx);
            end
        end
        met_time = toc;
        fprintf("finished excitation number %d, took %f minutes\n\n", excite, met_time/60)
    end
end
S = sum(S, 3);
t = 0:traj.dwellTime:traj.dwellTime*(size(traj.k_trajectory, 2)-1);
%apply T2 weighting
for excite = 1:size(S, 1)
    S(excite,:) = S(excite,:).*exp(-t/T2star);
end

S = MRSI_regrid(S, traj);
%convert to fid-a structure
out = MRSI_convert(S, traj, B0);

total = toc(func);
fprintf('total run time is %f minutes\n', total/60)
end

%FUNCTION readout: readout the signal from the spin system from all voxels
%get the readout from all voxels for metabolite m
function phantom_sig = readout(spins, Trc, scale, idx)
spins = spins(:, :, idx);
traced_spins = pagemtimes(spins, Trc);
sum_sig = sum(traced_spins, 3);
phantom_sig = scale*trace(sum_sig);
end


%FUNCTION calculate_H: create hamiltonians for each voxel based on gradient strength in
%grad matrix
function [H, H_inv] = calculate_H(phantom, gradients, m, time, B0, idx)
%get the new HAB for each spin at each voxel.
%This is derrived from ppm = ((v - v_ref)/v_ref)*10^6.
%
% (ppm*v_ref/10^6 + v_ref) = v. Where v_ref = (B0+G_x*x + G_y*y)*gamma
%
%(((ppm*v_ref/10^6 + v_ref)- v_ref_2)/v_ref_2)*10^6 = ((v-v_ref_2)/v_ref_2)*10^6 = dI_eff; where v_ref_2 = B0*gamma
%
%ppm_eff = (ppm*(B0+G*r)/10^6 + (B0+G*r) - B0)*10^6/B0
%(gamma can be factored out)


%get the metabolites
met = phantom.met(m);

%gyromagnetic ratio
gamma=42577000;  %[Hz/T]

gradients = gradients(:);

%compute the ppm_eff from above
ppm_eff = (gradients.*met.shifts(1)/1e6 + gradients - 1)*1e6;

%remove empty voxels
ppm_eff = ppm_eff(idx);
ppm_diff = ppm_eff - met.shifts(1);

%convert to rads/s
shift_rads = ppm_diff*2*pi*gamma*B0/1e6;
shift_rads = reshape(shift_rads, 1,1,length(shift_rads));

Fz = met.Fz;
Fz = repmat(Fz, [1,1,size(ppm_eff, 1)]); 


%element wise multiplication of shift rads. 
add_HAB = Fz.*shift_rads;

%sum along the third dimension (stack of matricies) and add HABJonly. 
%Final dimensions are (size of Hamiltonian along first dimension, size of
%Hamiltonian along second, and length of vectorized voxel positions)
HAB = met.HAB + add_HAB;
HAB = HAB*1i*time;
%initalize array for hamiltonian
H = zeros(size(Fz, 1), size(Fz, 2), sum(idx), 'single');

if(isdiag(HAB(:,:,1)))
    m = size(HAB,1);
    n = size(HAB,3);
    
    %get linearlized indexes of diagonals
    diag_idx = bsxfun(@plus, [1:m+1:m*m]', [0:n-1]*m*m);
    %extract diagonals from stack of matrices from idx
    diag = HAB(diag_idx);
    diag_exp = exp(diag);
    H(diag_idx) = diag_exp;
else
    for i = 1:size(H, 3)
        
        H(:,:, i) = expm(HAB(:,:,i));
    end
end
H_inv = conj(H);
end
