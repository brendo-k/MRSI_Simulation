
function [phantom] = MRSI_excite(phantom, deg, Fy, I0)
    angle = deg*(pi/180);
    HExcite = Fy*angle;
    for x = 1:size(phantom, 1)
        for y = 1:size(phantom, 2)
            if(~isequal(phantom(x,y).d, I0))
                phantom(x,y).d = expm(-1i*HExcite)*phantom(x,y).d*expm(1i*HExcite);
            end
        end
    end
end

