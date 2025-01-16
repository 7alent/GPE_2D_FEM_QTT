%
% This function is used in dmrg_eig(), it's written as an independent .m
% file for testing purpose
%

% Accumulates the right reduction W{k:d}'*A{k:d}*X{k:d}
function WAX1 = rightreduce_matrix(WAX2, w, A, x, rw1, n, rw2, Ra, ra1, ...
    ra2, rx1, m, rx2)
    % Right WAX has the form of the last matrix TT block, i.e. [ra, rw, rx]
    WAX1 = WAX2;
    wc = reshape(w, rw1, n*rw2);
    wc = conj(wc);
    xc = reshape(x, rx1*m, rx2);
    for k = 1:Ra
        WAX1{k} = reshape(WAX1{k}, ra2(k)*rw2, rx2);
        WAX1{k} = xc*WAX1{k}.'; % size rx1 m x ra2 rw2
        WAX1{k} = reshape(WAX1{k}, rx1, m*ra2(k)*rw2);
        WAX1{k} = WAX1{k}.';
        WAX1{k} = reshape(WAX1{k}, m*ra2(k), rw2*rx1);
        tmp = reshape(A{k}, ra1(k)*n, m*ra2(k));
        WAX1{k} = tmp*WAX1{k}; % size ra1(k)*n, rw2*rx1
        WAX1{k} = reshape(WAX1{k}, ra1(k), n*rw2*rx1);
        WAX1{k} = WAX1{k}.';
        WAX1{k} = reshape(WAX1{k}, n*rw2, rx1*ra1(k));
        WAX1{k} = wc*WAX1{k}; % size rw1, rx1 ra1
        WAX1{k} = reshape(WAX1{k}, rw1*rx1, ra1(k));
        WAX1{k} = WAX1{k}.';
    end
end