%
% This function is used in dmrg_eig(), it's written as an independent .m
% file for testing purpose
%

% A matrix-vectors product for the matrix in the 3D TT (WAX1-A-WAX2), and
% full vectors of size (rx1*m*rx2) x b. Returns (rw1*n*rw2) x b

%%% If WAX1 is 2-D rather than 3-D (e.g., when A is a unit matrix) this
%%% fucntion may not work due to wrong indexing

function w = local_matvec(x, rx1, m, rx2, b, rw1, n, rw2, WAX1, A, ...
    WAX2, Ra, ra1, ra2)
    xc = reshape(x, rx1*m*rx2, []);
    if isempty(b)
        b = size(xc, 2);
    end
    w = zeros(rw1*n*rw2, b);
    xc = xc.';
    xc = reshape(xc, b*rx1*m, rx2);
    for k = 1:Ra
        tmp = reshape(WAX2{k}, ra2(k)*rw2, rx2);
        wk = xc*tmp.';
        wk = reshape(wk, b*rx1, m*ra2(k)*rw2);
        wk = wk.';
        wk = reshape(wk, m*ra2(k), rw2*b*rx1);
        tmp = reshape(A{k}, ra1(k)*n, m*ra2(k));
        wk = tmp*wk;
        wk = reshape(wk, ra1(k)*n*rw2*b, rx1);
        wk = wk.';
        wk = reshape(wk, rx1*ra1(k), n*rw2*b);
        tmp = reshape(WAX1{k}, rw1, rx1*ra1(k));
        wk = tmp*wk;
        wk = reshape(wk, rw1*n*rw2, b);
        w = w+wk;
    end
end