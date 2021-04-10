clear 
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
%create phantom
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

excite = 1;
Kx = par.k_min_x:par.k_delta_x:par.k_max_x;
Ky = par.k_min_y:par.k_delta_y:par.k_max_y;
S = zeros(length(Kx), length(Ky), length(phan_x), length(phan_y), par.imageSize(3));


grad_t = zeros(size(gradient,2),1);
find_t = zeros(size(gradient,2),1);
indexing_t = zeros(size(gradient,2), 1);
itter_t = zeros(size(gradient,2), 1);
H_t = zeros(size(gradient,2)*64*64,1);
sand_t = zeros(size(gradient,2)*64*64, 1);
quer_t = zeros(size(gradient,2)*64*64, 1);


[grad_pos_x, grad_pos_y] = meshgrid(phan_x, phan_y); 
di_matrix = [phantom.dI];
di_matrix = reshape(di_matrix, [size(phantom, 1), size(phantom, 2)]);
trace_matrix = (Fx + 1i * Fy);
index_matrix = ones(size(phantom));

start=tic;
counter = 1;
for k=1:size(gradient,2)
    %Get current gradient and k space positions
    itteration_tic = tic;
    
    grad_tic = tic;
    cur_grad = gradient(excite,k).G;
    cur_time = gradient(excite,k).time;
    cur_k = k_traj(excite, k);
    grad_t(k) = toc(grad_tic);
    
    find_tic = tic;
    %Find k indecies for the readout point
    cur_Kx = find(Kx == real(cur_k), 1);
    cur_Ky = find(Ky == imag(cur_k), 1);
    find_t(k) = toc(find_tic);
    
    %Find the first index that is zero at Kx and Ky. This is where the singal
    %will be saved for thsi point
    index_tic = tic;
    spacial_index = index_matrix(cur_Kx, cur_Ky);
    indexing_t(k) = toc(index_tic);
    
    grad_matrix = real(cur_grad)*grad_pos_x + imag(cur_grad)*grad_pos_y + B0;
    freq_matrix = gamma*(di_matrix.*grad_matrix/1e6 + grad_matrix - B0)*2*pi;
    
    freq_map = containers.Map('KeyType','double','ValueType','any');

    for x = 1:size(phantom, 1)
        for y = 1:size(phantom, 2)
            
            
            if isKey(freq_map, freq_matrix(x,y))
                tic
                H = freq_map(freq_matrix(x,y));
                quer_t(counter) = toc;
            else
                tic
                times = (1i*freq_matrix(x,y)*cur_time);
                H = expm(times*Iz);
                freq_map(freq_matrix(x,y)) = H;
                H_t(counter) = toc;
            end
            
            temp= -1i*H*phantom(x,y).d*H;

            phantom(x,y).d = temp;
            S(cur_Kx, cur_Ky, x, y, spacial_index) = trace(trace_matrix * temp);
            counter = counter + 1;
        end
    end
    
    index_matrix(cur_Kx, cur_Ky) = index_matrix(cur_Kx, cur_Ky) + 1;
    itter_t(k) = toc(itteration_tic);
end
s_size = size(S);
S = reshape(S, [s_size(1), s_size(2), length(phan_x)*length(phan_y), par.imageSize(3)]);
S = sum(S, 3);
total_time = toc(start)

fprintf('grad %d\n', mean(grad_t))
fprintf('find %d\n', mean(find_t))
fprintf('indexing %d\n', mean(indexing_t))
fprintf('itteration %d\n', mean(itter_t))
fprintf('Hamil %d\n', mean(nonzeros(H_t)))
fprintf('query %d\n', mean(nonzeros(quer_t)))
fprintf('sandwich %d\n', mean(nonzeros(sand_t)))
