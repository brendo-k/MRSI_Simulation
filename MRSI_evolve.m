function [phantom] = MRSI_evolve(phantom, time, B0, Ix)
    gamma=42577000;  %[Hz/T]
    for x = 1:size(phantom, 1)
        for y = 1:size(phantom, 2)
            Hevol = (phantom(x,y).dI*B0*gamma*2*pi/10e6)*Ix;
            phantom(x,y).d = complex(expm(-1i*Hevol*time)*phantom(x,y).d*expm(1i*Hevol*time));
        end
    end
end

