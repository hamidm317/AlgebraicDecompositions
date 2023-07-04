function mat = rnd_Sym_matrix(n)

    tmp = rand(n, n, 1);
    
    mat = (tmp + tmp') / 2;

end

