
%Inputs: 
% Traj: trajectory class 
% gMax: Maximum gradient strength in mT/m
% B0: B0 magnetic field in T
function [out] = MRSI_simulate(traj, gMax, B0)

tic
if ~exist('B0', 'var')
    B0 = 3;
end

%TODO: Set up checmical shifts and coupling constants with sim_Hamiltonian
gamma=42577000;  %[Hz/T]
%TODO: need to use sim_Hamiltonian to use multiple spin systems
I0=complex([1 0;0 1]);
Ix=complex(0.5*[0 1;1 0]);
Iy=(1i/2)*[0 -1;1 0];
Iz=complex((1/2)*[1 0;0 -1]);
Fx = Ix;
Fy = Iy;

%Calculate gradient, k space, and spatical parameters
[gradient]  = load_trajectory(traj, gMax);

%Create phantom
phantom = MRSI_build_phantom(65, 65, 0.2, 0.2);

%create an matrix of x coordinates (used for speedup)
phan_x = [phantom.x];
phan_x = reshape(phan_x, [size(phantom, 1), size(phantom, 2)]);

%create an matrix of y coordinates (used for speedup)
phan_y = [phantom.y];
phan_y = reshape(phan_y, [size(phantom, 1), size(phantom, 2)]);

%create trace matrix for readout. (used for speedup)
trace_matrix = (Fx + 1i * Fy);

%create a matrix of the chemical shifts (used for speedups)
di_matrix = [phantom.dI];
di_matrix = reshape(di_matrix, [size(phantom, 1), size(phantom, 2)]);

%Initalize array for signal readout
S = zeros(size(traj.k_trajectory, 1), size(traj.k_trajectory,2));

S = complex(S, 0);

readout_length = size(gradient, 2);
x_phantom_size = size(phantom, 1);
y_phantom_size = size(phantom, 2);

%Now start readout:
phantom = MRSI_excite(phantom, 90, Fy, I0);
phantom = MRSI_evolve(phantom, 1e-3, 3, Ix);
phantom = MRSI_excite(phantom, 180, Fy, I0);
parfor excite=1:size(gradient,1)
    fprintf("simulating excitation number %d\n", excite)
    
    %Excite phantom
    new_phantom = phantom;
    
    for k=1:readout_length
        %Get current gradient and k space positions
        cur_grad = gradient(excite,k).G;
        cur_time = gradient(excite,k).time;
      
        grad_matrix = real(cur_grad)*phan_x + imag(cur_grad)*phan_y + B0;
        
        %calcuate the frequency of each position
        freq_matrix = gamma*(di_matrix.*grad_matrix/1e6 + grad_matrix - B0)*2*pi;
        
        phantom_sig = zeros(size(new_phantom));
        %loop through phantom
        for x = 1:x_phantom_size
            for y = 1:y_phantom_size
                
                %if nothing is there skip
                if isequal(new_phantom(x,y).d, I0)
                    continue;
                else
                    %get the number of rotations
                    rotations = (freq_matrix(x,y)*cur_time);
                    
                    %create the hamiltonian
                    matrix_exp = expm(-1i*rotations*Iz);
                    %inverse hamiltonian is just the complex conjugate
                    inverse_exp = expm(1i*rotations*Iz);
                    
                    %apply sandwich operation
                    new_phantom(x,y).d = matrix_exp*new_phantom(x,y).d*inverse_exp;
                    
                    if k~=1
                        %save to signal
                        phantom_sig(x,y) = trace(trace_matrix * new_phantom(x,y).d);
                    end
                end
            end
        end
        
        if k == 1
            new_phantom = MRSI_evolve(new_phantom, 1e-3-cur_time, 3, Ix);
            for x = 1:x_phantom_size
                for y = 1:y_phantom_size
                    
                    %if nothing is there skip
                    if isequal(new_phantom(x,y).d, I0)
                        continue;
                    else
                        phantom_sig(x,y) = trace(trace_matrix * new_phantom(x,y).d);
                    end
                end
            end
        end
        
        S(excite, k) = sum(phantom_sig, 'all');
        %index_matrix(cur_Kx, cur_Ky) = spacial_index + 1;
    end
end
S = MRSI_regrid(S, traj);
%TODO: apply t2 weighting for entire signal

%convert to fid-a structure
out = MRSI_convert(S, traj, B0);
toc
end

