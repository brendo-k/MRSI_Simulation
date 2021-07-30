function [phantom] = MRSI_excite(phantom, deg, direction)
    %direction should be either x  or y
    if(~strcmp(direction, 'y') && ~strcmp(direction, 'x'))
        error('direciton should be either x or y');
    end

    %convert flip angle to readians
    flip_angle = deg*(pi/180);
    for x = 1:size(phantom, 1)
        for y = 1:size(phantom, 2)
            for m = 1:length(phantom(x,y).met)
                %getting proper excitation hamiltonian for flip direction
                if(strcmp(direction, 'y'))
                    HExcite = phantom(x, y).met(m).Fy*flip_angle;
                else
                    HExcite = phantom(x, y).met(m).Fx*flip_angle;
                end
                %sandwich opperator
                phantom(x,y).d{m} = expm(-1i*HExcite)*phantom(x,y).d{m}*...
                    expm(1i*HExcite);
            end
        end
    end
end


