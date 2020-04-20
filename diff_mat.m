function D = diff_mat(m, order)
% dimension : m
% order     : order (default value is 2)

    if order == 4
        D = diag(sparse(2/3*ones(1, m-1)), 1) + diag(sparse(-2/3*ones(1, m-1)),-1) ...
            + diag(sparse(-1/12*ones(1, m-2)),2) + diag(sparse(1/12*ones(1, m-2)),-2);
        D(1, end-1:end) = [1/12, -2/3]; D(2, end) = 1/12; 
        D(end-1:end, 1) = [-1/12, 2/3]; D(end, 2) = -1/12;
    else
        D = diag(sparse(1/2*ones(1, m-1)),1) + diag(sparse(-1/2*ones(1, m-1)),-1);
        D(1, m) = -1/2; D(m, 1) = 1/2;
    end
end