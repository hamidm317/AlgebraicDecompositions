function res = diag_zeroer(A)

    n = min(size(A));
    res = A;

    for i = 1 : n

        res(i, i) = 0;

    end

end

