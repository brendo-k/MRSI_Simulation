function b0Slice = MRSI_getB0SliceFromNii(fileName, sliceNumber, deltaT)
    arguments
        fileName (1, :) char {mustBeFile}
        sliceNumber (1, 1) double
        deltaT (1, 1) double
    end
    niftiStruct = spm_vol(fileName);
    [b0Map, ~] = spm_read_vols(niftiStruct);
    b0Map = flip(permute(b0Map, [2,1,3]), 1);
    deltaPhaseSlice = b0Map(:, :, sliceNumber);
    b0Slice = (deltaPhaseSlice*(pi/180))/(getGamma('overTwoPi', false) * deltaT);
end