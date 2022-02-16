function AB = mmx_mult(A,B,mod)
%MMX_MULT It generalizes "mmx" by Yuval on MATLAB File Exchange to complex
%numbers.
%
% =========================================================================
% Author: Yi-Hao Chen
% email: jollycyclopes@gmail.com
% date: 6/6/2018
%
% Acknowledgement:
%   This code is based on Yural's "mmx".
%

% The calculation supports only "double" precision.
Atype = class(A);
input_A_is_single = false;
if isequal(Atype,'gpuArray')
    Atype = classUnderlying(A);
end
if isequal(Atype,'single')
    input_A_is_single = true;
    A = double(A);
end
Btype = class(B);
input_B_is_single = false;
if isequal(Btype,'gpuArray')
    Btype = classUnderlying(B);
end
if isequal(Btype,'single')
    input_B_is_single = true;
    B = double(B);
end

if isreal(A) && isreal(B)
    if exist('mod','var')
        AB = mmx('mult',A,B,mod);
    else
        AB = mmx('mult',A,B);
    end
else
    realA = real(A);
    imagA = imag(A);

    realB = real(B);
    imagB = imag(B);

    if exist('mod','var')
        AB = mmx('mult',realA,realB,mod) - mmx('mult',imagA,imagB,mod)...
         + 1i*( mmx('mult',realA,imagB,mod) + mmx('mult',imagA,realB,mod) );
    else
        AB = mmx('mult',realA,realB) - mmx('mult',imagA,imagB)...
         + 1i*( mmx('mult',realA,imagB) + mmx('mult',imagA,realB) );
    end
end

% Transform it back to "single" precision if using single precision.
if input_A_is_single || input_B_is_single
    AB = single(AB);
end

end