
function simulate_MRSI(k_traj, dwell_time, gMax)
    close all;
    

    %TODO: Set up checmical shifts and coupling constants with sim_Hamiltonian
    dI=-3.0; %water chemical shift[ppm]
    
    %Calculate gradient, k space, and spatical parameters
    [k_traj, gradient, k_par] = load_trajectory(k_traj, dwell_time, gMax);
    
    k_max_x = max(real(k_traj), 'all');
    k_max_y = max(imag(k_traj), 'all');
    k_min_x = min(real(k_traj), 'all');
    k_min_y = min(imag(k_traj), 'all');

    k_fov_x = k_max_x - k_min_x;
    k_fov_y = k_max_y - k_min_y;

    delta_x = 1/k_fov_x;
    delta_y = 1/k_fov_y;

    %Other parameters;
    np=64; %[spacial points]
    tp=64; %[spectral points]
    
    fov = 0.2; %[m]
    deltaX = fov/np; %[m]
    B0=3; %[T]
    gradient=20e-3; %[mT/m]
    gamma=42577000;  %[Hz/T]
    sw = fov*gradient*gamma; %[Hz] Derrivation: fov = 1/delta_k = 1/(gradient * gamma * delta_t) -> fov*gradient*gamma = 1/delta_t = sw
    dwelltime=1/sw; %[s] Derrivation: delta_t = 1/sw
    t=[0:dwelltime:dwelltime*(np-1)]; %[s]
    deltaK = 1/fov; %[m^-1]
    fovk = 1/deltaX; %[m^-1]
    sw_spectral=1/(t(end)*2); %[Hz]
    dt_spectral=1/sw_spectral; %[Hz] dwell time of the spectral dimension
    t_spectral=[0:dt_spectral:dt_spectral*(tp-1)]; %[s]
    
    %BK No echo time needed
    %te=0.030; %30 ms echo time;
    
    %set up the basis matrices (Calling the two spins 'I' and 'S', but they
    %could also be called 'A' and 'B' or whatever):
    %first the 1-spin basis:
    I0=[1 0;0 1];
    Ix=0.5*[0 1;1 0];
    Iy=(1i/2)*[0 -1;1 0];
    Iz=(1/2)*[1 0;0 -1];
    
    %S0=I0;
    %Sx=Ix;
    %Sy=Iy;
    %Sz=Iz;
    
    %BK - Don't need 2 spin basis
    % I0Sx=kron(I0,Sx);
    % I0Sy=kron(I0,Sy);
    % I0Sz=kron(I0,Sz);
    % 
    % IxS0=kron(Ix,S0);
    % IyS0=kron(Iy,S0);
    % IzS0=kron(Iz,S0);
    
    
    %Now set up the F matrices:
    % Fx=IxS0+I0Sx;
    % Fy=IyS0+I0Sy;
    % Fz=IzS0+I0Sz;
    
    %Since only one spin basis
    Fx = Ix;
    Fy = Iy;
    
    %set up the IS matrix:
    %IS=IxS0*I0Sx + IyS0*I0Sy + IzS0*I0Sz;
    
    %Now set up the equilibrium Density matrix:
    d0=Iz;
    
    %Set up the hamiltonian operators:
    H90=Iy*pi/2;
    %H180=Ix*pi;
    
    %Set up the free evolution hamiltonian operators for the delay periods:
    %Hevol_dS=(dS*B0*gamma*2*pi/1e6)*I0Sz;
    %Hevol_JIS=(JIS*2*pi)*IS;
    %Hevol=Hevol_dI+Hevol_dS+Hevol_JIS;
    
    %Now do the pulse sequence:
    dtemp=expm(-1i*H90)*d0*expm(1i*H90);                 %90 degree excitation
    
    %phatom parameters
    phantom_size = 0.2;
    phantom_deltaX = 0.2/64;
    
    %x coordinates for the phantom
    x = -phantom_size/2 + phantom_deltaX/2:phantom_deltaX:phantom_size/2 - phantom_deltaX/2;
    
    %create phantom
    for i = 1:64
    
        %Set coordinate
        phantom(i).x = x(i);
    
        %water from 20 to 40
        if i >= 20 && i <= 40 
            phantom(i).d = dtemp;
            phantom(i).dI = dI;
        else
            %Nothing elsewhere
            phantom(i).d = I0;
            phantom(i).dI = 0;
        end
    end
    
    %calculate k coordinates
    Kx = -fovk/2+deltaK/2:deltaK:fovk/2-deltaX/2;
    k_start = Kx(1);
    
    %Calculate time it takes to get to first k space position
    %Uses equation k=gamma*G*t which becomes t = k/(gamma*G)
    evolution_time = k_start/(gamma*-gradient);
    
    %Move to first k space position (gradient echo)
    for i = 1:64
        %Get gradient for position x
        b_grad = -gradient*phantom(i).x;
    
        %calculate the frequency at position x(j)
        freq_encode = phantom(i).dI*(B0+b_grad)*gamma/1e6 + gamma*(B0+b_grad) - gamma*B0; %(v-vref)*1e6/vref
    
        %Get evolution hamiltonian
        Hevol_start=(freq_encode*2*pi)*Iz;
    
        %Sandwhich opperation for the duration of evolution_time
        phantom(i).d = expm(-1i*Hevol_start*evolution_time)*phantom(i).d*expm(1i*Hevol_start*evolution_time);
    end
    
    %aquire signal from readout
    S = zeros(length(t), tp);
    
    %Now start readout:
    for k=1:tp
        for n=1:length(t)
            for j = 1:length(phantom)
                %get signal
                S(n,k) = S(n,k) + trace((Fx + 1i * Fy) * phantom(j).d);
    
                %calculate the gradient at position x(j)
                b_grad = gradient*phantom(j).x;
    
                %calculate the frequency at position x(j)
                freq_encode = phantom(j).dI*(B0+b_grad)*gamma/1e6 + gamma*(B0+b_grad) - gamma*B0;
    
                %hamiltonian and evolve
                Hevol1=(freq_encode*2*pi)*Iz;
                phantom(j).d = expm(-1i*Hevol1*dwelltime)*phantom(j).d*expm(1i*Hevol1*dwelltime);
            end
        end
    
        %rewind gradients
        for j = 1:length(phantom)
    
            %calculate the gradient at position x(j)
            b_grad = -gradient*phantom(j).x;
    
            %calculate the frequency at position x(j)
            freq_encode = phantom(j).dI*(B0+b_grad)*gamma/1e6 + gamma*(B0+b_grad) - gamma*B0;
    
            %hamiltonian and evolve
            Hevol1=(freq_encode*2*pi)*Iz;
            phantom(j).d = expm(-1i*Hevol1*(t(end)+dwelltime))*phantom(j).d*expm(1i*Hevol1*(t(end)+dwelltime));
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
    
    
    
    %PART  B - TAKE THE FFT AND PLOT THE SPECTRUM:
    
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
