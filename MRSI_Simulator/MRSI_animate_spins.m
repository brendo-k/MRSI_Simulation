function MRSI_animate_spins(voxels)

    figure
    for excite = 1:size(voxels, 1)        
        [x,y] = meshgrid(1:size(voxels,3), 1:size(voxels,4));
        ramping = ramping_interpolate(squeeze(voxels(excite, 1, :, :)), 20);
        vox_angles = squeeze(ramping(1, :, :));
        u = real(vox_angles);
        v = imag(vox_angles);
        h = quiver(x, y, u, v, 0.25);
        for t = 1:size(ramping, 1)
            vox_angles = squeeze(ramping(t, :, :));
            h.UData = real(vox_angles);
            h.VData = imag(vox_angles);
            pause(0.05)
        end
        xout = interp1(1:size(voxels,2), squeeze(voxels(excite, :,:,:)), 1:0.1:size(voxels,2), 'spline'); 
        for t = 1:100
            vox_angles = squeeze(xout(t, :, :));
            h.UData = real(vox_angles);
            h.VData = imag(vox_angles);
            pause(0.05);
        end
    end
end

function points = ramping_interpolate(first, n)
    start_pos = complex(ones(size(first, 1), size(first,2)), 0);
    start_pos(first == 0) = 0;
    start_angles = angle(start_pos);
    end_angles = angle(first);
    diff_angle = end_angles - start_angles;
    delta_angle = diff_angle/(n-1);
    delta_angle = repmat(delta_angle, [1,1,n]);
    multiply = 0:n-1;
    for i = 1:n
        delta_angle(:,:,i) = delta_angle(:,:,i).*multiply(i);
    end
    points = repmat(start_pos, [1,1,n]);
    points = points .* exp(1i.*delta_angle);
    points = permute(points, [3,1,2]);
end
