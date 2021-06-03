
%Inputs:
% Traj: trajectory class
% gMax: Maximum gradient strength in mT/m
% B0: B0 magnetic field in T
% T2star: T2star weighting in seconds
function [out, voxel_sig] = MRSI_simulate(traj, phantom, gMax, B0, T2star)

tic
if ~exist('phantom', 'var')
    disp('hello')
    metabolites = cell(64,64);
    metabolites(10:50, 10:50) = sysLac;
    metabolites(20:40, 20:40) = sysH2O;
    %Create phantom
    phantom = MRSI_build_phantom([0.2, 0.2], metabolites, B0);
end
if ~exist('B0', 'var')
    B0 = 3;
end
if ~exist('gMax', 'var')
    gMax = 10;
end

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
phantom = MRSI_excite(phantom, 90, 'y');
phantom = MRSI_evolve(phantom, 1e-3);
phantom = MRSI_excite(phantom, 180, 'x');

voxel_sig = complex(zeros(size(gradient,1), size(gradient, 2), x_phantom_size, y_phantom_size), 0);

parfor excite=1:size(gradient,1)
    fprintf("simulating excitation number %d\n", excite)
    
    %Excite phantom
    new_phantom = phantom;
    
    for k=1:readout_length
        %Get current gradient and k space positions
        cur_grad = gradient(excite,k).G;
        cur_time = gradient(excite,k).time;
        
        grad_matrix = real(cur_grad)*phan_x + imag(cur_grad)*phan_y + B0;
        phantom_sig = zeros(size(new_phantom));
        %loop through phantom
        for x = 1:x_phantom_size
            for y = 1:y_phantom_size
                for m = 1:length(new_phantom(x,y).met)
                    dI_eff = (new_phantom(x,y).met(m).shifts*grad_matrix(x,y)/1e6 + grad_matrix(x,y) - B0)*1e6/B0;
                    shift_rads = dI_eff*-2*pi*gamma*B0/1e6;
                    Iz = new_phantom(x,y).met(m).Iz ;
                    HAB = complex(zeros(size(Iz,1), size(Iz,2)),0);
                    for spin  = 1:size(Iz,3)
                        HAB = HAB + Iz(:,:,spin)*shift_rads(spin);
                    end
                    HAB = HAB + new_phantom(x,y).met(m).HABJonly;
                    %create the hamiltonian
                    matrix_exp = expm(-1i*HAB*cur_time);
                    %inverse hamiltonian is just the complex conjugate
                    inverse_exp = expm(1i*HAB*cur_time);
                    
                    %apply sandwich operation
                    new_phantom(x,y).met(m).d{1} = matrix_exp*new_phantom(x,y).met(m).d{1}...
                        *inverse_exp;
                   
                    if k~=1
                        %save to signal
                        trace_matrix = new_phantom(x,y).met(m).Fx + 1i*new_phantom(x,y).met(m).Fy;
                        phantom_sig(x,y) = trace(trace_matrix * new_phantom(x,y).met(m).d{1});
                        voxel_sig(excite, k, x, y) = phantom_sig(x,y);
                    end
                end
                
            end
        end
        
        if k == 1
            new_phantom = MRSI_evolve(new_phantom, 1e-3-cur_time);
            for x = 1:x_phantom_size
                for y = 1:y_phantom_size
                    for m = 1:length(new_phantom(x,y).met)
                        trace_matrix = new_phantom(x,y).met(m).Fx + 1i*new_phantom(x,y).met(m).Fy;
                        phantom_sig(x,y) = trace(trace_matrix * new_phantom(x,y).met(m).d{1});
                        voxel_sig(excite, k, x, y) = phantom_sig(x,y);
                    end
                end
            end
        end
        
        S(excite, k) = sum(phantom_sig, 'all');
    end
    
end

t = 0:1/traj.sw:1/traj.sw*(traj.imageSize(3)-1);
%apply T2 weighting
for excite = 1:size(S, 1)
    S(excite,:) = S(excite,:).*exp(-t/T2star);
end

S = MRSI_regrid(S, traj);
%TODO: apply t2 weighting for entire signal

%convert to fid-a structure
out = MRSI_convert(S, traj, B0);

toc
end

