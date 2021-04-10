
function [specs, S] = simulate_MRSI(k_traj, par, gMax)
    close all;

    %TODO: Set up checmical shifts and coupling constants with sim_Hamiltonian
    dI=-3.0; %water chemical shift[ppm]
    
    %Calculate gradient, k space, and spatical parameters
    [k_traj, gradient, par] = load_trajectory(k_traj, par, gMax);
    
    %Other parameters;
    %np=64; %[spacial points]
    %tp=64; %[spectral points]
    
    fovX = par.fovX; %[m]
    deltaX = par.delta_x; %[m]
    fovY = par.fovY;
    deltaY = par.delta_y;
    X = -fovX/2 + deltaX/2:deltaX:fovX/2 - deltaX/2;
    Y = -fovY/2 + deltaY/2:deltaY:fovY/2 - deltaY/2;
    
    B0=3; %[T] 
    gamma=42577000;  %[Hz/T]
    
    %set up the basis matrices (Calling the two spins 'I' and 'S', but they
    %could also be called 'A' and 'B' or whatever):
    %first the 1-spin basis:
    %TODO: need to use sim_Hamiltonian to use multiple spin systems
    I0=[1 0;0 1];
    Ix=0.5*[0 1;1 0];
    Iy=(1i/2)*[0 -1;1 0];
    Iz=(1/2)*[1 0;0 -1];

    Fx = Ix;
    Fy = Iy;
    
    %Now set up the equilibrium Density matrix:
    d0=Iz;
    
    %phatom parameters
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
    for i = 1:length(phan_x)
        for j = 1:length(phan_y)
            %Set coordinate
            phantom(i,j).x = phan_x(i);
            phantom(i,j).y = phan_y(j);

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

    
    %calculate k coordinates
    k_fov_x = par.k_fov_x;
    k_fov_y = par.k_fov_y;
    k_delta_x = par.k_delta_x;
    k_delta_y = par.k_delta_y;
    Kx = par.k_min_x:k_delta_x:par.k_max_x;
    Ky = par.k_min_y:k_delta_y:par.k_max_y;
    
    %aquire signal from readout
    S = zeros(length(Kx), length(Ky), par.imageSize(3));
    
    %Now start readout:
    for excite=1:size(gradient,1)
        excite_tic = tic
        %Reset phantom to transverse magnetization
        phantom = MRSI_reset(phantom, d0);
        %Excite phantom
        phantom = MRSI_excite(phantom, 90);

        for k=1:size(gradient,2)
            %Get current gradient and k space positions
            cur_grad = gradient(excite,k).G;
            cur_time = gradient(excite,k).time;
            cur_k = k_traj(excite, k);

            %Find k indecies for the readout point
            cur_Kx = find((Kx==real(cur_k)),1);
            cur_Ky = find((Ky==imag(cur_k)),1);

            %Find the first index that is zero at Kx and Ky. This is where the singal 
            %will be saved for thsi point
            spacial_index = find(S(cur_Kx, cur_Ky, :));
            phantom_time = tic;
            for x = 1:size(phantom,1)
                for y = 1:size(phantom,2)
                    %Get the magentic field x and y coponents at position x and y 
                    grad_x = real(cur_grad)*phantom(x,y).x;
                    grad_y = imag(cur_grad)*phantom(x,y).y;
                    
                    tic
                    %calculate the frequency at position x(j)
                    freq_encode = phantom(x,y).dI*(B0+grad_x+grad_y)*gamma/1e6 + gamma*(B0+grad_x+grad_y) - gamma*B0;
                    freq_cal = toc;
                    
                    tic
                    %hamiltonian and evolve
                    Hevol1=(freq_encode*2*pi)*Iz;
                    phantom(x,y).d = expm(-1i*Hevol1*cur_time)*phantom(x,y).d*expm(1i*Hevol1*cur_time);
                    Hamiltonian = toc;
                    
                    tic
                    signal = trace((Fx + 1i * Fy) * phantom(x,y).d);
                    trace_t = toc;

                    tic
                    %get signal
                    S(cur_Kx, cur_Ky,spacial_index) = S(cur_Kx, cur_Ky, spacial_index) + signal;
                    readout = toc;
                end
            end
            t_phantom = toc(phantom_time);
            %disp('time to calculate phantom: ' + string(t_phantom));
            
        end
        toc(excite_tic)
    end
    
    %Never applied T2 weighting 
    %TODO: apply t2 weighting for entire signal
    %T2star=0.1; %[s]
    %y = exp(-t/T2star);  %apply T2* weighting;
    %S = S.*y'
    
    specs=fftshift(fft(S,1),1);
    specs=fftshift(fft(S,2),2);
    specs=fftshift(fft(S,3),3);
    
%     %Set up spectral frequency axis (done the same way as spacial frequency)
%     f_spec= -par.sw/2:par.sw/(par.imageSize(3)-1):par.sw/2;
%     
%     %x coord calculated with Fov = 1/delta_k 
%     %                            = 1/(gamma*G*delta_t)
%     %                            = SW/(gamma*G)    
%     %Convert frequency axis to ppm:
%     ppm=f_spec/(B0*gamma/1e6);
%     
%     
%     %Plot the spectrum:
%     figure;
%     if(size(spec, 2) == 1)
%         plot(x_coord, abs(spec))
%     else
%         surf(ppm, x_coord, abs(spec(:,20,:)), 'LineWidth',1.2);
%     end
%     ylabel('x coordinates','FontSize',20);
%     xlabel('ppm', 'Fontsize', 20)
%     zlabel('Spectral Intensity (a.u.)','FontSize',20);
%     title('Solution Figure for Part B','FontSize',24);
%     box off;
%     set(gcf,'color','w');

end

