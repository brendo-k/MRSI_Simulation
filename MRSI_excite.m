function phantom = MRSI_excite(phantom, deg)
    rad = deg*pi/180;
    Iy=(1i/2)*[0 -1;1 0];
    H_excite=Iy*rad;

    for x = 1:size(phantom, 1)
        for y = 1:size(phantom, 2)
            phantom(x,y).d = expm(1i*H_excite)*phantom(x,y).d*expm(1i*H_excite);
        end
    end
end
