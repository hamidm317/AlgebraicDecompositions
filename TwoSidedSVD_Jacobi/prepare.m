function [U_, V_, A_] = prepare(A)

    [m, n] = size(A);
    
    if m < n
        
        A = A';
        
        [V_, A_] = qr(A);
        U_ = eye(m);
        A_ = A_';
        
    else
        
        [U_, A_] = qr(A);
        V_ = eye(n);
        
        
        
    end

end

