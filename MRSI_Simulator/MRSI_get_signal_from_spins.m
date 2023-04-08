%FUNCTION readout: readout the signal from the spin system from all voxels
%get the readout from all voxels for metabolite m
function scaled_signal = MRSI_get_signal_from_spins(phantom)
    scaled_signal = 0;
    for m = 1:length(phantom.met)
        scale = 2 ^ (2 - phantom.met(m).nspins);
        trc = single(phantom.met(m).Fx + 1i*phantom.met(m).Fy);
        sumSpins = sum(phantom.spins{m}, [3, 4]);
        sumSpins = squeeze(sumSpins);
        signal = trace(sumSpins*trc);
        scaled_signal = scaled_signal + scale*signal;
    end
end