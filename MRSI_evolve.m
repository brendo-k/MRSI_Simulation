function [phantom] = MRSI_evolve(phantom, time, B0)
    gamma=42577000;  %[Hz/T]
    for x = 1:size(phantom, 1)
        for y = 1:size(phantom, 2)
            for m = 1:length(phantom(i,j).met)
                Ix = phantom(x,y).met(m).Fx;
                Hevol = (phantom(x,y).met(m).HAB*B0*gamma*2*pi/10e6)*Ix;
                phantom(x,y).met(m).d = expm(-1i*Hevol*time)*phantom{x,y}.d*...
                                        expm(1i*Hevol*time);
            end
        end
    end
end

