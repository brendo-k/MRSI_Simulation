function S_remap = MRSI_regrid(S, traj)
    coord1 = cell2mat(traj.K_coordinates.k(1));
    coord2 = cell2mat(traj.K_coordinates.k(2));
    S_remap = zeros(length(coord1), length(coord2), size(S,2));
    k_traj = traj.k_trajectory;
    
    for readout = 1:size(S,1)
        if(length(unique(k_traj(readout,:))) == 1)
            %cartesian MRSI
            k_y = imag(k_traj(readout,1));
            k_x = real(k_traj(readout,1));
            
            x_i = (coord1 == k_x);
            y_i = (coord2 == k_y);
            S_remap(x_i, y_i, :) = S(readout,:);
        end
    end
    
end
