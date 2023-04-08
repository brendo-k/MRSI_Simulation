phantom = MRSI_MNI_phantom(3, 50);

trajectory = phaseEncoded(imageSize=[8, 8]);
signal = MRSI_simulate_gpu(trajectory, phantom);
save('simulated_signal', signal)
ft_final = op_CSIFourierTransform(signal);
save('ft_final_signal', ft_final)
op_CSIPlot(ft_final)