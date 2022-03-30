function MRSI_make_sparse(phantom, threshold)
    arguments
        phantom (1, 1) struct
        threshold (1, 1) double = 1e-5
    end
    
    for i = numel(phantom.spins)
        spinMatrix = phantom.spins{i};
        spinMatrix(spinMatrix < threshold) = 0;
        phantom.spins{i} = spinMatrix;
    end
end