function [phantom] = MRSI_excite(phantom, deg, direction)
    %direction should be either x  or y
    if(~strcmp(direction, 'y') && ~strcmp(direction, 'x'))
        error('direciton should be either x or y');
    end
    
    %convert flip angle to readians
    flip_angle = deg*(pi/180);
    
    for m = 1:length(phantom.met)
        for y = 1:length(phantom.y)
            for x = 1:length(phantom.x)
                %getting proper excitation hamiltonian for flip direction
                if(strcmp(direction, 'y'))
                    HExcite = phantom.met(m).Fy*flip_angle;
                else
                    HExcite = phantom.met(m).Fx*flip_angle;
                end
                %sandwich opperator
                phantom.spins{m}(y,x,:,:) = expm(-1i*HExcite)*squeeze(phantom.spins{m}(y,x,:,:))*...
                    expm(1i*HExcite);
            end
        end
    end
    
end


