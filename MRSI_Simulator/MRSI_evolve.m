%MRSI_evolve.m
%Lets the spins of the phantom go thorugh free evolution for a given time
%Input
%   phantom:Phantom object from MRSI_build_phantom
%   time:   time fo evolution (s)
%   B0:     Magnetic field strength b0;
%
%Output
%   phantom: phantom after free evolution 

function [phantom] = MRSI_evolve(phantom, time)
    gamma=42577000;  %[Hz/T]
    for x = 1:size(phantom, 1)
        for y = 1:size(phantom, 2)
            for m = 1:length(phantom(x,y).met)
                Hevol = phantom(x,y).met(m).HAB;
                phantom(x,y).d{m} = expm(-1i*Hevol*time)*phantom(x,y).d{m}*...
                                        expm(1i*Hevol*time);
            end
        end
    end
end

