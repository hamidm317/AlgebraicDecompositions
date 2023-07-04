function L = lambda_generator(lambda, P)

    L = zeros(P, P, P);
    
    for i = 1 : P
        
        L(i, i, i) = lambda(i);
        
    end

end

