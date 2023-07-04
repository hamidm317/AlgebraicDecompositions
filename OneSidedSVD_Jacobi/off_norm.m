function off = off_norm(A)

    tmp_ = diag_zeroer(A);
    off = norm(tmp_, 'fro');

end

