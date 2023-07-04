function [G, U] = HOOI_Tucker(T, R)

    G = [];
    U = [];
    
    if length(size(T)) == 3 && length(R) == 3
        
        U0 = U_initializer(T, R);
        
%         U1 = cell2mat(U0{1});
        U2 = cell2mat(U0(2));
        U3 = cell2mat(U0(3));
        
        do_it = 1;
        max_iter = 100;
        iteration = 0;
        while do_it
            
%             Z_temp = tmprod(T, {U2', U3'}, [2, 3]);
            Z_temp = ttm(T, {U2', U3'}, [2, 3]).data;
            U1 = Un_step(Z_temp, R(1), 1);
            
%             Z_temp = tmprod(T, {U1', U3'}, [1, 3]);
            Z_temp = ttm(T, {U1', U3'}, [1, 3]).data;
            U2 = Un_step(Z_temp, R(2), 2);
            
%             Z_temp = tmprod(T, {U1', U2'}, [1, 2]);
            Z_temp = ttm(T, {U1', U2'}, [1, 2]).data;
            U3 = Un_step(Z_temp, R(3), 3);
            
            if iteration > max_iter
                
                do_it = 0;
                
            end
            
            iteration = iteration + 1;
            
        end
        
%         G = tmprod(T, {U1', U2', U3'}, 1 : 3);
        G = ttm(T, {U1', U2', U3'}, 1 : 3).data;
        U = {U1, U2, U3};
        

        
    end

end

