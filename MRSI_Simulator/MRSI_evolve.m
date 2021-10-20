%MRSI_evolve.m
%Lets the spins of the phantom go thorugh free evolution for a given time
%Input
%   phantom:Phantom object from MRSI_build_phantom
%   time:   time fo evolution (s)
%   B0:     Magnetic field strength b0;
%   MemoryOptions:
%       'use_disc', 1 or 0
%   ArgumentOptions:
%       'argument_type', 'struct' or 'matrix'
%
%Output
%   phantom: phantom after free evolution 

function [phantom] = MRSI_evolve(phantom, time, MemoryOptions, ArgumentOptions)
arguments
    phantom 
    time (1,1) double {mustBeNonnegative}
    MemoryOptions.use_disc (1,1) logical = 0
    ArgumentOptions.argument_type {mustBeMember(ArgumentOptions.argument_type,{'struct', 'matrix'})} = "struct"
    ArgumentOptions.HAB (:,:) double
end

tic
%STUCT PASSED INTO PHANTOM VARIABLE
if(strcmp(ArgumentOptions.argument_type, 'struct'))
    for m = 1:length(phantom.met)
        if(MemoryOptions.use_disc)
            spins = load_spins(phantom, m);
        else
            spins = phantom.spins{m};
        end
        Hevol = expm(phantom.met(m).HAB*time*-1i);
        inv_Hevol = expm(phantom.met(m).HAB*time*1i);

        spins = calculate(spins, Hevol, inv_Hevol);

        if(MemoryOptions.use_disc)
            phantom.file{m} = save_spins(spins, phantom.met_names{m}, 'MRSI_evolve');
        else
            phantom.spins{m} = spins;
        end
        clear spins;
    end
%MATRIX PASSED INTO STRUCTURE VARIABLE
elseif(strcmp(ArgumentOptions.argument_type, 'matrix'))
    Hevol = expm(ArgumentOptions.HAB*time*-1i);
    inv_Hevol = expm(ArgumentOptions.HAB*time*1i);
    phantom = calculate(phantom, Hevol, inv_Hevol);
end
cpu = toc;
fprintf('total time from cpu: %6.2f\n', cpu)
end

function spins = calculate(spins, H, H_inv)
    phan_size = size(spins, [1,2]);
    spins = permute(spins, [3,4,1,2]);
    spins = reshape(spins, size(spins,1), size(spins,2), []);
    spins = pagemtimes(pagemtimes(H_inv,spins), H);
    spins = reshape(spins, size(spins,1), size(spins,2), phan_size(1), phan_size(2));
    spins = permute(spins, [3,4,1,2]);
end

% tic
% for m = 1:length(phantom.met)
%     if(MemoryOptions.use_disc)
%         spins = load_spins(phantom.file{m}, phantom, m);
%     else
%         spins = phantom.spins{m};
%     end    
%     Hevol = expm(phantom.met(m).HAB*time*-1i);
%     inv_Hevol = expm(phantom.met(m).HAB*time*1i);
%     spins = gpuArray(spins);
%     spins = permute(spins, [3,4,1,2]);
%     spins = reshape(spins, size(spins,1), size(spins,2), []);
%     spins = pagemtimes(pagemtimes(Hevol,squeeze(spins)), inv_Hevol);
%     spins = reshape(spins, size(spins,1), size(spins,2), length(phantom.y), length(phantom.x));
%     spins = permute(spins, [3,4,1,2]);
%     spins = gather(spins);
%     
%     if(MemoryOptions.use_disc)
%         phantom.file{m} = save_spins(spins, phantom.met_names{m}, 'MRSI_evolve');
%     else
%         phantom.spins{m} = spins;
%     end
%     clear spins;
% end
% toc
% gpu = toc;
% fprintf('total time from gpu: %f\n', gpu)

