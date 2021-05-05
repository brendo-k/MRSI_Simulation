function [out] = MRSI_zero_fill(in, n)
    zero_fill = zeros(in.sz(1), in.sz(2) ,n);
    out = in;
    out.fids = cat(3, in.fids, zero_fill);
    
end