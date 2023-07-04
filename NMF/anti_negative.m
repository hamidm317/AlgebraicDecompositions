function nn_B = anti_negative(B)

    [m, n] = size(B);
    
    nn_B = B;
    
    for i = 1 : m * n
        
        if B(i) < 0
            
            nn_B(i) = 0;
            
        end
        
    end

end

