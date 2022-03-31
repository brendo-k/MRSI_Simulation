%MRSI_evolve.m
%Lets the spins of the phantom go thorugh free evolution for a given time
%Input
%   phantom:Phantom object from MRSI_build_phantom
%   time:   time fo evolution (s)
%   B0:     Magnetic field strength b0;
%   MemoryOptions:
%       'use_disc', 1 or 0
%
%Output
%   phantom: phantom after free evolution

function [phantom] = MRSI_excite(phantom, deg, direction, MemoryOptions, ArgumentOptions)
    arguments
        phantom
        deg (1,1) double
        direction (:,1) char
        MemoryOptions.use_disc (1,1) logical = 0
        ArgumentOptions.argument_type {mustBeMember( ...
            ArgumentOptions.argument_type, {'struct', 'matrix'})} = "struct"
        ArgumentOptions.F (:,:) double
    end
    if(strcmp(ArgumentOptions.argument_type, 'struct'))
        %direction should be either x  or y
        if(~strcmp(direction, 'y') && ~strcmp(direction, 'x'))
            error('direciton should be either x or y');
        end

        %convert flip angle to readians
        flip_angle = deg*(pi/180);

        for m = 1:length(phantom.met)
            if(MemoryOptions.use_disc)
                spins = load_spins(phantom, m);
            else
                spins = phantom.spins{m};
            end
            %getting proper excitation hamiltonian for flip direction
            if(strcmp(direction, 'y'))
                HExcite = expm(1i*phantom.met(m).Fy*flip_angle);
                inv_HExcite = expm(-1i*phantom.met(m).Fy*flip_angle);
            else
                HExcite = expm(1i*phantom.met(m).Fx*flip_angle);
                inv_HExcite = expm(-1i*phantom.met(m).Fx*flip_angle);
            end

            spins = calculate(spins, HExcite, inv_HExcite);
            %sandwich operation

            if(MemoryOptions.use_disc)
                phantom.file{m} = save_spins(spins, phantom.met_names{m}, 'MRSI_excite');
            else
                phantom.spins{m} = spins;
            end

        end
    elseif(strcmp(ArgumentOptions.argument_type, 'matrix'))
        flip_angle = deg*(pi/180);
        HExcite = expm(-1i*ArgumentOptions.F*flip_angle);
        inv_HExcite = expm(1i*ArgumentOptions.F*flip_angle);
        phantom = calculate(phantom, HExcite, inv_HExcite);
    end
end

function spins = calculate(spins, H, H_inv)
    imageSize = size(spins, [3, 4]);
    spins = reshape(spins, size(spins, 1), size(spins, 2), []);
    spins = pagemtimes(pagemtimes(H, spins), H_inv);
    spins = reshape(spins, [size(spins, [1, 2]) imageSize]);
end


