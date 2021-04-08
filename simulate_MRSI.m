
function simulate_MRSI(k_traj, dwell_time, gMax)
    close all;
    

    %TODO: Set up checmical shifts and coupling constants with sim_Hamiltonian
    dI=-3.0; %water chemical shift[ppm]
    
    %Calculate gradient, k space, and spatical parameters
    [k_traj, gradient, scan_par] = load_trajectory(k_traj, dwell_time, gMax);
    
    %Other parameters;
    %np=64; %[spacial points]
    %tp=64; %[spectral points]
    
    fovX = scan_par.fovX; %[m]
    deltaX = fovX.delta_x; %[m]
    fovY = scan_par.fovY;
    deltaY = scan_par.delta_y;
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
    
    %Set up the hamiltonian operators:
    H90=Iy*pi/2;
    
    %Excite spin density
    dtemp=expm(-1i*H90)*d0*expm(1i*H90);                 %90 degree excitation
    
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
    
    %calculate k coordinates
    k_fov_x = par.k_fov_x;
    k_fov_y = par.k_fov_y;
    k_delta_x = par.k_delta_x;
    k_delta_y = par.k_delta_y;
    Kx = -k_fov_x/2+k_delta_x/2:k_delta_x:k_fov_x/2-k_delta_x/2;
    Ky = -k_fov_y/2+k_delta_y/2:k_delta_y:k_fov_y/2-k_delta_y/2;
    
    %aquire signal from readout
    S = zeros(size(k_traj));
    
    %Now start readout:
    for excite=1:size(gradient,1)
        
        for k=1:size(gradient,2)
            for i = 1:size(phantom,1)
                for j = 1:size(phantom,2)
                    %calculate the gradient at position x(j)
                    cur_grad = gradient(excite,k).G;
                    grad_x = real(cur_grad)*phantom(i,j).x;
                    grad_y = imag(cur_grad)*phantom(i,j).y;
                    
                    %calculate the frequency at position x(j)
                    freq_encode = phantom(j).dI*(B0+grad_x+grad_y)*gamma/1e6 + gamma*(B0+grad_x+grad_y) - gamma*B0;
                    
                    %hamiltonian and evolve
                    Hevol1=(freq_encode*2*pi)*Iz;
                    phantom(j).d = expm(-1i*Hevol1*dwelltime)*phantom(j).d*expm(1i*Hevol1*dwelltime);
                    
                    %get signal
                    S(excite,k) = S(n,k) + trace((Fx + 1i * Fy) * phantom(j).d);
                end
            end
        end
    end
    
    %Never applied T2 weighting 
    %TODO: apply t2 weighting for entire signal
    %T2star=0.1; %[s]
    %y = exp(-t/T2star);  %apply T2* weighting;
    %S = S.*y'
    
    %plot the FID:
    hold on;
    for i = 1:size(S,2)
        %plotting the fid of each k space acquisition
        plot(t,S(:,i),'LineWidth',1);
    end
    xlabel('Time (s)','FontSize',20);
    ylabel('FID signal intensity (a.u.)','FontSize',20);
    title('Solution Figure for Part A','FontSize',24);
    box off;
    hold off;
    set(gcf,'color','w');
    
    %zero padding
    %S= vertcat(S, zeros(64,1));
    %Apply 2FFT to FID in order to get spectral and spacial dimensions
    spec=fftshift(fft2(S));
    
    %Set up spacial frequency axis
    f=[(-sw/2)+(dwelltime/2):1/((t(end)-t(1))+dwelltime):(sw/2)-(dwelltime/2)];
    
    %Set up spectral frequency axis (done the same way as spacial frequency)
    f_spec=[(-sw_spectral/2)+(dt_spectral/2):1/((t_spectral(end)-t_spectral(1))+dt_spectral):(sw_spectral/2)-(dt_spectral/2)];
    
    %x coord calculated with Fov = 1/delta_k 
    %                            = 1/(gamma*G*delta_t)
    %                            = SW/(gamma*G)
    x_coord=f/(gamma*gradient);
    
    %Convert frequency axis to ppm:
    ppm=f_spec/(B0*gamma/1e6);
    
    
    %Plot the spectrum:
    figure;
    if(size(spec, 2) == 1)
        plot(x_coord, abs(spec))
    else
        surf(ppm, x_coord, abs(spec), 'LineWidth',1.2);
    end
    ylabel('x coordinates','FontSize',20);
    xlabel('ppm', 'Fontsize', 20)
    zlabel('Spectral Intensity (a.u.)','FontSize',20);
    title('Solution Figure for Part B','FontSize',24);
    box off;
    set(gcf,'color','w');

end
