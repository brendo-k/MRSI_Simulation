function phantom = MRSI_reset(phantom, d0, I0)
    for x = 1:size(phantom, 1)
        for y = 1:size(phantom, 2)
            %check if something is there and reset it to Iz state
            if(~isequal(phantom(x,y).d, I0))
                phantom(x,y).d = d0;
            end
        end
    end
end
