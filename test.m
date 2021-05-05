

tic
B0 = 3;
gMax = 3;


%TODO: Set up checmical shifts and coupling constants with sim_Hamiltonian
gamma=42577000;  %[Hz/T]
%TODO: need to use sim_Hamiltonian to use multiple spin systems
I0=complex([1 0;0 1]);
Ix=complex(0.5*[0 1;1 0]);
Iy=(1i/2)*[0 -1;1 0];
Iz=complex((1/2)*[1 0;0 -1]);
Fx = Ix;
Fy = Iy;

%Now set up the equilibrium Density matrix:
d0=Iz;

%Calculate gradient, k space, and spatical parameters
[gradient]  = load_trajectory(traj, gMax);

%Create phantom
phantom = MRSI_build_phantom(64, 64, 0.2, 0.2);

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

k_first = traj.K_coordinates.k{1};
k_second = traj.K_coordinates.k{2};

%Initalize array for signal readout
S = zeros(length(k_first), length(k_second), size(phantom, 1), ... 
    size(phantom, 2), traj.imageSize(3));

S = complex(S, 0);

k_traj = traj.k_trajectory;

%Now start readout:
for excite=1:size(gradient,1)
    fprintf("simulating excitation number %d\n", excite)
    
    %Reset phantom to transverse magnetization
    phantom = MRSI_reset(phantom, d0, I0);
    %Excite phantom
    phantom = MRSI_excite(phantom, 90, Fy, I0);

    %this matrix gives us the index   spectral point to be saved for readout 
    index_matrix = ones(size(phantom));
  
    for k=1:size(gradient,2)
        %Get current gradient and k space positions
        cur_grad = gradient(excite,k).G;
        cur_time = gradient(excite,k).time;
        cur_k = k_traj(excite, k);
        
        %Find k indecies for the readout point
        cur_Kx = find((k_first==real(cur_k)),1);
        cur_Ky = find((k_second==imag(cur_k)),1);
        
        %Find the first index that is zero at Kx and Ky. This is where the singal
        %will be saved for thsi point
        spacial_index = index_matrix(cur_Kx, cur_Ky);
        
        %calculate a matrix of the magentic field at each i,j position of
        %the phantom
        grad_matrix = real(cur_grad)*phan_x + imag(cur_grad)*phan_y + B0;
        
        %calcuate the frequency of each position
        freq_matrix = gamma*(di_matrix.*grad_matrix/1e6 + grad_matrix - B0)*2*pi;
        
        %loop through phantom
        for x = 1:size(phantom, 1)
            for y = 1:size(phantom, 2)
                
                %if nothing is there skip
                if isequal(phantom(x,y).d, I0)
                    continue;
                else
                    %get the number of rotations
                    rotations = (freq_matrix(x,y)*cur_time);
                    
                    %create the hamiltonian
                    matrix_exp = expm(-1i*rotations*Iz);
                    %inverse hamiltonian is just the complex conjugate
                    inverse_exp = conj(matrix_exp);
                    
                    %apply sandwich operation
                    phantom(x,y).d = matrix_exp*phantom(x,y).d*inverse_exp;
                    %save to signal
                    S(cur_Kx, cur_Ky, x, y, spacial_index) = trace(trace_matrix * phantom(x,y).d);
                end
            end
        end
        
        index_matrix(cur_Kx, cur_Ky) = spacial_index + 1;
    end
end

s_size = size(S);
%reshape to create one dimension for x and y dimension of phantom
S = reshape(S, [s_size(1), s_size(2), length(phan_x)*length(phan_y), s_size(5)]);
%sum the signal from x and y dimensions of the phantom
S = squeeze(sum(S, 3));

%TODO: apply t2 weighting for entire signal

%convert to fid-a structure
out = MRSI_convert(S, traj, B0);
toc
end

