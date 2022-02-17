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
function [out, spin_animation] = MRSI_simulate_gpu(traj, phantom, scanParameters, MemoryOptions, Debug)
    arguments
        traj (1,1) Trajectory
        phantom (:, :) struct
        scanParameters.gMax (1,1) double = 30
        scanParameters.B0 (1,1) double = 3
        scanParameters.B0_map (:, :) double = zeros(length(phantom.y), length(phantom.x));
        MemoryOptions.use_disc (1,1) logical = 0;
        Debug.plot_spins (1, 1) logical = 0;
        Debug.spinEcho (1, 1) logical = false;
        Debug.plot_spins_tr (1, 1) double = 0;
    end
    functionTimer = tic;
    gMax = scanParameters.gMax;
    B0 = scanParameters.B0;
    B0_map = scanParameters.B0_map;

    %if the phantom is empty, data is stored in files.
    if(isempty(phantom.spins{1}))
        MemoryOptions.use_disc = 1;
    end
   
    %Calculate gradients and time steps for each gradient
    [gradient, gradientTime] = MRSI_load_ktrajectory(traj, gMax);
    gradientTime = gradientTime;

    %Initalize array for signal readout. Dims are (num Metabolites, TR, readout points)
    MRSISignal = zeros([length(phantom.met), size(gradient, [1, 2])]);
    MRSISignal = complex(MRSISignal, 0);

    isCartesian = false;
    if(strcmp(traj.name, 'cartesian')); isCartesian = true; end
    TE = 0.010;
    %voxel_sig = complex(zeros(size(gradient,1), size(gradient, 2), x_phantom_size, y_phantom_size), 0);
    if(Debug.plot_spins == true)
        %TODO: save spins to visualize
    end
    for m = 1:length(phantom.met)
        spins = getSpins(MemoryOptions, phantom, m);
        spins = MRSI_excite(spins, 90, 'y', 'argument_type', 'matrix', 'F', phantom.met(m).Fy);
        if(Debug.spinEcho)
            spins = MRSI_evolve(spins, TE, 'argument_type', 'matrix', 'HAB', phantom.met(m).HAB);
            spins = MRSI_excite(spins, 180, 'x', 'argument_type', 'matrix', 'F', phantom.met(m).Fx);
        end
        [spins, nonZeroIndex] = vectorizeSpins(spins);
        fprintf("simulating metabolite %s\n", phantom.met_names{m})

        met = phantom.met(m);
        %change variables to be gpuArrays
        [trc, Iz, HAB] = convertToGPU(met);
        [xGrid, yGrid] = getVoxelCoordinateGrid(phantom);

        metaboliteSignal = complex(zeros(size(gradient, 1, 2)), 0);
        shieldingFactor = convertToShieldingFactor(met.shifts);
        for trNumber = 1:size(gradient, 1)
            trTimer = tic;
            fprintf("simulating TR number %d\n", trNumber)

            %set spins for this TR.
            trSpins = spins;
            trSpins = gpuArray(trSpins);
            timeElapsed = 0;
            for k=1:size(gradient, 2)
                %get map of gradient for each voxel position in the grid
                curGradient = gradient(trNumber, k);
                gradientMap = getGradientMap(curGradient, xGrid, yGrid, B0_map, B0);
                gradientMap = gradientMap(squeeze(nonZeroIndex));
                %Calculate hamiltonians
                timeStep = gradientTime(trNumber, k);
                timeElapsed = timeElapsed + timeStep;
                if(~isCartesian) || k <= 2
                    gradientFrequency = getFrequencyFromGradients(gradientMap, shieldingFactor, met, B0);
                    HAB_effective = calculateNewHAB(Iz, gradientFrequency, HAB);

                    [H, H_inv, is_diag] = calculateHamiltonians(HAB_effective, timeStep);
                end
                %Apply Hamiltonians
                trSpins = applyHamiltonians(trSpins, H, H_inv, is_diag);
                if(k == 1 && Debug.spinEcho)
                    timeElapsed = timeElapsed - timeStep;
                    trSpins = MRSI_evolve(trSpins, TE - timeStep, 'argument_type', 'matrix', 'HAB', phantom.met(m).HAB);
                end
                %get signal
                scale = 2 ^ (2 - met.nspins);
                metaboliteSignal(trNumber, k) = getSignalFromSpins(trSpins, trc, scale);
                
                if(Debug.plot_spins == true && Debug.plot_spins_tr == trNumber)
                    
                end
                
            end
            trExecutionTime = toc(trTimer);
            fprintf("finished TR number %d, took %f minutes\n", ...
                trNumber, trExecutionTime/60)
        end
        MRSISignal(m, :, :) = applyT2Decay(traj, metaboliteSignal, phantom, m);
    end
    %add up all the metabolite signals.
    MRSISignal = sum(MRSISignal, 1);
    %Change the rotation direction
    MRSISignal = conj(MRSISignal);

    MRSISignal = MRSI_regrid(MRSISignal, traj);
    %convert to fid-a structure
    out = MRSI_convert(MRSISignal, traj, B0);

    total = toc(functionTimer);
    fprintf('total run time is %f minutes\n', total/60)
end


%FUNCTION readout: readout the signal from the spin system from all voxels
%get the readout from all voxels for metabolite m
function scaledSignal = getSignalFromSpins(spins, Trc, scale)
    
    sumSpins = sum(spins, 3);
    signal = trace(sumSpins*Trc);
    scaledSignal = scale*signal;
end


%FUNCTION calculate_H: create hamiltonians for each voxel based on gradient strength in
function [H, H_inv, is_diag] = calculateHamiltonians(new_HAB, time)
    %get the new HAB for each spin at each voxel.

    if(isdiag(new_HAB(:,:,1)))
        m = size(new_HAB, 1);
        n = size(new_HAB, 3);
        diag_idx = get_diag_index(m,n);

        HAB_evolve = new_HAB(diag_idx).*(1i*time);
        %initalize array for hamiltonian

        %get linearlized indexes of diagonals
        %extract diagonals from stack of matrices from idx
        exp_HAB = exp(HAB_evolve);
        exp_conj = conj(exp_HAB);
        H = exp_HAB;
        H_inv = exp_conj;
        is_diag = true;
    else
        HAB_evolve = new_HAB.*(1i*time);
        try
            H = myexpm_(HAB_evolve, [], [], false, false, true);
        catch exception
            disp('Spins system too big! Splitting into two')
            if (strcmp(exception.identifier,'parallel:gpu:array:OOM'))
                half = floor(size(HAB_evolve, 3)/2);
                H1 = myexpm_(HAB_evolve(:, :, 1:half), [], [], false, false, true);
                H2 = myexpm_(HAB_evolve(:, :, half+1:end), [], [], false, false, true);
                H = cat(3, H1, H2);
            else
                rethrow(exception);
            end
        end
        H_inv = conj(H);
        is_diag = false;
    end
end



function spins = getSpins(MemoryOptions, phantom, m)
    if(MemoryOptions.use_disc)
        spins = load_spins(phantom, m);
    else
        spins = phantom.spins{m};
    end
end



function [spins, nonZeroIndex] = vectorizeSpins(spins)
    %vectorize the position dimension. Now excite spins is a stack of
    %matricies
    spins = reshape(spins, size(spins,1), size(spins,2), []);
    %excite_spins = gpuArray(excite_spins);
    nonZeroIndex = any(spins, [1,2]);
    spins = spins(:,:, nonZeroIndex);
end



function gradientMap = getGradientMap(gradient, xGrid, yGrid, B0Map, b0)
    %apply gradient to x and y directions to form a gradient matrix.
    %grad_matrix has gradient strength at the x and y position;
    currentGradientX = real(gradient);
    currentGradientY = imag(gradient);
    gradientMap = currentGradientX*xGrid + currentGradientY*yGrid;
    gradientMap = gradientMap + B0Map;
end

function shieldingFactor = convertToShieldingFactor(shifts)
    shieldingFactor = (shifts)/1e6;
end



function spins = applyHamiltonians(spins, H, Hinv, isDiag)
    %if it is diagonal we can multiply matrix by a vector of
    %diagonals instead of a matrix with diagonals.
    if(isDiag)
        HReshaped = reshape(H, 1, size(H,1), size(H,2));
        HInvReshaped = reshape(Hinv, size(H,1), 1, size(H,2));
        spins = spins.*HReshaped;
        spins = spins.*HInvReshaped;
    else
        spins = pagemtimes(pagemtimes(Hinv, spins), H);
    end
end



function signalWithDecay = applyT2Decay(traj, metaboliteSignal, phantom, m)
    t = traj.t;
    signalWithDecay = metaboliteSignal .* exp(-t/phantom.T2(m));
end


function [xGrid, yGrid] = getVoxelCoordinateGrid(phantom)
    %create an matrix of x and y coordinates (used for speedup)
    [xGrid, yGrid] = meshgrid(single(phantom.x), single(phantom.y));
end

function [trc,  Iz, HAB] = convertToGPU(met)
    trc = gpuArray(single(met.Fx + 1i*met.Fy));
    Iz = gpuArray(single(met.Iz));
    HAB = gpuArray(single(met.HAB));
end

function offsetFrequency = getFrequencyFromGradients(gradients, shielding, met, b0)
    %gyromagnetic ratio
    gamma = -getGamma('overTwoPi', false);
    ppm = met.shifts;
    offsetFrequency = gamma * gradients * (ppm' / 1e6 + 1);
%     currentFrequency = shielding * gradients * gamma;
%     offsetFrequency = currentFrequency - gamma * b0 - met.shifts_rads;
end

%calculate new HAB. We can take advantage that HAB is already calculated. We
%only need to add Iz_1*shift1Gradient + IZ_2*shift2Gradient
function new_HAB = calculateNewHAB(Iz, gradientFrequency, HAB)
    Iz = permute(Iz, [1,3,2]);
    IzMultiplied = pagefun(@mtimes, Iz, gradientFrequency');
    IzMultiplied = permute(IzMultiplied, [1,3,2]);
    
    %sum along the third dimension (stack of matricies) and add HABJonly.
    %Final dimensions are (size of Hamiltonian along first dimension, size of
    %Hamiltonian along second, and length of vectorized voxel positions)
    new_HAB = repmat(HAB, [1, 1, size(IzMultiplied, 3)]);
    new_HAB = new_HAB + IzMultiplied;
end