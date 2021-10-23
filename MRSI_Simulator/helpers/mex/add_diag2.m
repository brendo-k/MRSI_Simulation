function H = add_diag2(H, v)
    diag_idx = get_diag_index(size(H, 1), size(H,3));
    H(diag_idx) = H(diag_idx) + v;
end