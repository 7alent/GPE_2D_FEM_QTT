%
% This function is used in dmrg_eig(), it's written as an independent .m
% file for testing purpose
%

% Accumulates the left reduction W{1:k}'*A{1:k}*X{1:k}
function WAX2 = leftreduce_matrix(WAX1, w, A, x, rw1, n, rw2, Ra, ra1, ...
    ra2, rx1, m, rx2)
    % Left WAX has the form of the first matrix TT block, i.e. [rw, rx, ra]
    WAX2 = WAX1;
    wc = reshape(w, rw1, n*rw2);
    xc = reshape(x, rx1*m, rx2);
    for k = 1:Ra
        WAX2{k} = reshape(WAX2{k}, rw1, rx1*ra1(k));
        WAX2{k} = wc'*WAX2{k}; % size n rw2 x rx1 ra1
        WAX2{k} = reshape(WAX2{k}, n, rw2*rx1*ra1(k));
        WAX2{k} = WAX2{k}.';
        WAX2{k} = reshape(WAX2{k}, rw2*rx1, ra1(k)*n);
        tmp = reshape(A{k}, ra1(k)*n, m*ra2(k));
        WAX2{k} = WAX2{k}*tmp; % size rw2 rx1 m ra2
        WAX2{k} = reshape(WAX2{k}, rw2, rx1*m*ra2(k));
        WAX2{k} = WAX2{k}.';
        WAX2{k} = reshape(WAX2{k}, rx1*m, ra2(k)*rw2);
        WAX2{k} = xc.'*WAX2{k}; % size rx2, ra2 rw2
        WAX2{k} = reshape(WAX2{k}, rx2*ra2(k), rw2);
        WAX2{k} = WAX2{k}.';
    end
end