%% start
clear 
start = tic;
gradient = load('grad.mat').grad;
k_traj = load('traj.mat').traj;
par = load('par.mat').par;

I0=[1 0;0 1];
Ix=0.5*[0 1;1 0];
Iy=(1i/2)*[0 -1;1 0];
Iz=(1/2)*[1 0;0 -1];
d0 = Iz;
dI=-3.0;
B0 = 3;
gamma=42577000;

Fx = Ix;
Fy = Iy;

phantom_nx = 64;
phantom_ny = 64;
phantom_sizeX = 0.2; %[m]
phantom_sizeY = 0.2; %[m]
phantom_deltaX = phantom_sizeX/phantom_nx;
phantom_deltaY = phantom_sizeY/phantom_ny;

%x coordinates for the phantom
phan_x = -phantom_sizeX/2 + phantom_deltaX/2:phantom_deltaX:phantom_sizeX/2 - phantom_deltaX/2;
phan_y = -phantom_sizeY/2 + phantom_deltaY/2:phantom_deltaY:phantom_sizeY/2 - phantom_deltaY/2;

phantom = struct('d', cell(length(phan_x), length(phan_y)), 'dI', cell(length(phan_x),length(phan_y)));
%% create phantom
for i = 1:length(phan_y)
    for j = 1:length(phan_x)
        %Set coordinate
        phantom(i,j).x = phan_x(j);
        phantom(i,j).y = phan_y(i);
        
        %water from 20 to 40, already excited
        if i >= 20 && i <= 44 && j >= 20 && j<=44
            phantom(i,j).d = d0;
            phantom(i,j).dI = dI;
        else
            %Nothing elsewhere
            phantom(i,j).d = I0;
            phantom(i,j).dI = 0;
        end
    end
end

%% initalize variables
excite = 1;
Kx = par.k_min_x:par.k_delta_x:par.k_max_x;
Ky = par.k_min_y:par.k_delta_y:par.k_max_y;
S = zeros(length(Kx), length(Ky), length(phan_x), length(phan_y), par.imageSize(3));

[grad_pos_x, grad_pos_y] = meshgrid(phan_x, phan_y); 
di_matrix = [phantom.dI];
di_matrix = reshape(di_matrix, [size(phantom, 1), size(phantom, 2)]);
trace_matrix = (Fx + 1i * Fy);
index_matrix = ones(size(phantom));
a1 = zeros(size(gradient,2)*400, 1);
a2 = zeros(size(gradient,2)*400, 1);
counter = 1;
%% scan
for k=1:size(gradient,2)
    %Get current gradient and k space positions
    cur_grad = gradient(excite,k).G;
    cur_time = gradient(excite,k).time;
    cur_k = k_traj(excite, k);
    
    %Find k indecies for the readout point
    cur_Kx = find(Kx == real(cur_k), 1);
    cur_Ky = find(Ky == imag(cur_k), 1);
    
    %Find the first index that is zero at Kx and Ky. This is where the singal
    %will be saved for thsi point
    spacial_index = index_matrix(cur_Kx, cur_Ky);
    
    grad_matrix = real(cur_grad)*grad_pos_x + imag(cur_grad)*grad_pos_y + B0;
    freq_matrix = gamma*(di_matrix.*grad_matrix/1e6 + grad_matrix - B0)*2*pi;
    
    %freq_map = containers.Map('KeyType','double','ValueType','any');
    for x = 1:size(phantom, 1)
        for y = 1:ceil(size(phantom, 2)/2)
            if isequal(phantom(x,y).d, I0)
                continue;
            else
                
                times = (freq_matrix(x,y)*cur_time);
                
                matrix_exp = expm(-1i*times*Iz);
                inverse_exp = conj(matrix_exp);
                
                phantom(x,y).d = matrix_exp*phantom(x,y).d*inverse_exp;
                S(cur_Kx, cur_Ky, x, y, spacial_index) = trace(trace_matrix * phantom(x,y).d);

                if(y ~= ceil(size(phantom, 2)/2) || mod(size(phantom,2),1) ~= 1)
                    new_y = size(phantom,2)+1-y;
                    phantom(x,new_y).d = inverse_exp*phantom(x,new_y).d*matrix_exp;
                    S(cur_Kx, cur_Ky, x, new_y, spacial_index) = trace(trace_matrix * phantom(x,new_y).d);
                end
                counter = counter + 1;
                
            end
        end
    end
    
    index_matrix(cur_Kx, cur_Ky) = index_matrix(cur_Kx, cur_Ky) + 1;
end
s_size = size(S);
S = reshape(S, [s_size(1), s_size(2), length(phan_x)*length(phan_y), par.imageSize(3)]);
S = sum(S, 3);
toc(start);
