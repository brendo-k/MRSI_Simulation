function X = multslash(M,L)
%MULTSLASH Compute L(:,:,...)/M(:,:,...) for each 2D slice of an array (M) 
%          with arbitrary dimensions support.
%
% INPUT
%   M  : n_D array (m1 x n x [p x q x ...]), for all possible m1,n=1,2,3,... 
%        and optional higher dimensions.
%   L  : n_D array (m2 x n x [p x q x ...]), for all possible m2,n=1,2,3,... 
%        and optional higher dimensions.
%   
% OUTPUT
%   X  : n_D array (m2 x m1 x [p x q x  ...]).
%
% NOTE 1 -- This function may use a large amount of memory for huge array. 
%           Test before usage.
%
% NOTE 2 -- Underdetermined system (more unknowns than equation)
%   The solution is basic solution obtained with sparse mldivide
%   which is not the same as basic solution when calling for full matrix.
%
% See also: multiprod
%
% Author: Xiaodong Qi <i2000s@hotmail.com>
% History: original 26-Apr-2011
% Inspired by Bruno Luong's MultiSolver code to solve a linear equ system.
%
% Copyright (c) 2011, Xiaodong Qi
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
% =========================================================================
% Author: Yi-Hao Chen
% email: jollycyclopes@gmail.com
% date: 6/6/2018
%
% Acknowledgement:
%   This code is based on "multinv" by Xiaodong Qi.
%

% The calculation supports only "double" precision.
Mtype = class(M);
input_M_is_single = false;
if isequal(Mtype,'gpuArray')
    Mtype = classUnderlying(M);
end
if isequal(Mtype,'single')
    input_M_is_single = true;
    M = double(M);
end
Ltype = class(L);
input_L_is_single = false;
if isequal(Ltype,'gpuArray')
    Ltype = classUnderlying(L);
end
if isequal(Ltype,'single')
    input_L_is_single = true;
    L = double(L);
end

sM = size(M);
sL = size(L);
[m1,n1] = deal(sM(1),sM(2));
[m2,n2] = deal(sL(1),sL(2));
if n1 ~= n2
   error('multslash: Their column dimension should be the same.');
else
    n = n1;
end
if ~isequal(sM(3:end),sL(3:end))
    error('multslash: The dimensions after the 2nd one should be the same for both input matrices.');
end
nM = prod(sM(3:end));
M = reshape(M,[m1,n,nM]);

% Build sparse matrix and solve
MI = reshape(1:m1*nM,m1,1,nM);
MI = repmat(MI,[1 n 1]); % m1 x n x nM
MJ = reshape(1:n*nM,1,n,nM);
MJ = repmat(MJ,[m1 1 1]); % m1 x n x nM

M = sparse(MI(:),MJ(:),M(:));
L = cell2mat(permute(mat2cell(reshape(L,m2,n,nM),m2,n,ones(nM,1)),[1 3 2]));

clear MI MJ
X = L / M; % The result will be different from L/full(M) if M isn't composed of square matices, but it's still correct.
clear L M
X = reshape(X,[m2,m1,sM(3:end)]);

% Transform it back to "single" precision if using single precision.
if input_M_is_single || input_L_is_single
    X = single(X);
end

end
