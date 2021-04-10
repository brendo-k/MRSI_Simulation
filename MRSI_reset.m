function phantom = MRSI_reset(phantom, d0)
    I0=[1 0;0 1];
    for x = 1:size(phantom, 1)
        for y = 1:size(phantom, 2)
            %If nothing is there
            if(~isequal(phantom(x,y).d, I0))
                phantom(x,y).d = d0;
            end
        end
    end
end
