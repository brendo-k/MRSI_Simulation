function [spinSignal, spinCoordinate] = MRSI_show_spins(spin, xGrid, yGrid, met)
    vectX = xGrid(:);
    vectY = yGrid(:);
    spinCoordinate = [vectX, vectY];
    traceMatrix = (met.Fx + 1i*met.Fy);
    traceSpin = pagemtimes(spin,traceMatrix);
    spinSignal = zeros(size(vectX));
    for iMatrix = 1:prod(size(traceSpin, [3,4]))
            spinSignal(iMatrix) = trace(squeeze(traceSpin(:, :, iMatrix)));
    end
    quiver(spinCoordinate(:,1), spinCoordinate(:,2), real(spinSignal), imag(spinSignal))
end
