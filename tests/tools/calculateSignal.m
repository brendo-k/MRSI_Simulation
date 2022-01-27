
function signal = calculateSignal(phantom, spins)
    Fx = phantom.Fx;
    Fy = phantom.Fy;
    signal = trace(pagemtimes((Fx + 1i*Fy), spins));
end