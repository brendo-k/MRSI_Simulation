function S_remap = MRSI_regrid(S, traj)
arguments
    S (:,:) double
    traj (1,1) Trajectory
end


if(contains(traj.name, 'cart', 'IgnoreCase', 1))
    %cartesian MRSI
    S_remap = zeros(traj.imageSize);
    S_remap = permute(S_remap, [3,1,2]);
    
    for readout = 1:size(S,1)
        [x, y] = ind2sub(traj.imageSize(1:2),readout);
        S_remap(:, x, y) = S(readout,:);
    end
    S_remap = permute(S_remap, [1,3,2]);
else
    S_remap = permute(S, [2, 1]);
end

    
end
