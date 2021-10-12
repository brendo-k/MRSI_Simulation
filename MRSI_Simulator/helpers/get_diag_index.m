% get_diag_index.m
%
% Returns linearlized indexes of the diagonal in a page of matricies
%
%
%

function idx = get_diag_index(dim1_len, dim3_len)
    numel_matrix = dim1_len*dim1_len;
    %gives 
    page_idx = 1:dim1_len+1:numel_matrix;
    page_offset = (0:dim3_len-1)*numel_matrix;
    idx = page_idx' + page_offset;
end