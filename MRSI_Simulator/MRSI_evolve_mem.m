%MRSI_evolve.m
%Lets the spins of the phantom go thorugh free evolution for a given time
%Input
%   phantom:Phantom object from MRSI_build_phantom
%   time:   time fo evolution (s)
%   B0:     Magnetic field strength b0;
%
%Output
%   phantom: phantom after free evolution 

function [phantom] = MRSI_evolve_mem(phantom, time)
tic
for m = 1:length(phantom.met)
    file_id = fopen(phantom.file{m});
    spins = fread(file_id,inf, 'single=>single');
    fclose(file_id);
    
    spins = reshape(spins, [length(phantom.y), length(phantom.x), size(phantom.d{m}, [1,2])]);
    Hevol = expm(phantom.met(m).HAB*time*-1i);
    inv_Hevol = expm(phantom.met(m).HAB*time*1i);
    
    spins = permute(spins, [3,4,1,2]);
    spins = reshape(spins, size(spins,1), size(spins,2), []);
    spins = pagemtimes(pagemtimes(Hevol,squeeze(spins)), inv_Hevol);
    spins = reshape(spins, size(spins,1), size(spins,2), length(phantom.y), length(phantom.x));
    spins = permute(spins, [3,4,1,2]);
    phantom.spins{m} = spins;
    clear spins;
end
toc
end

