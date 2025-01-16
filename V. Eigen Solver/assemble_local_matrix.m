%
% This function is used in dmrg2(), it's written as an independent .m
% file for testing purpose
%

% Builds the full (rw1*n*rw2) x (rx1*m*rx2) matrix from its TT blocks
function [B, sparseflag] = assemble_local_matrix(WAX1, A, WAX2, Ra, ...
    ra1, ra2, rw1, n, rw2, rx1, m, rx2)
    % Check the sparsity of the matrix blocks
    sparseflag = true;
    for k = 1:Ra
        if ~issparse(A{k})
            sparseflag = false;
        end
    end
    if sparseflag
        B = sparse(rw2*rw1*n, rx2*rx1*m); % reverse order !!!
        % The reverse order is needed since usually the n x m part is large
        % and sparse, so let it be the senior dimension.
        % Note that currently only canonical sparse matrices are allowed
        for k = 1:Ra
            tmp = reshape(WAX2{k}, rw2, rx2);
            tmp = sparse(tmp);
            Bk = reshape(WAX1{k}, rw1, rx1);
            Bk = sparse(Bk);
            Bk = kron(Bk, tmp); % mind endiannes
            Bk = kron(A{k}, Bk); % mind endiannes
            B = B+Bk;
        end
    else
        % There are dense blocks, everything is dense, and in the natural 
        % index order
        B = zeros(rw1*n*rw2, rx1*m*rx2);
        for k = 1:Ra
            Bk = reshape(WAX1{k}, rw1*rx1, ra1(k));
            tmp = reshape(A{k}, ra1(k), n*m*ra2(k));
            if issparse(tmp)
                % Don't mix sparse if we are already full
                tmp = full(tmp);
            end
            Bk = Bk*tmp;
            Bk = reshape(Bk, rw1, rx1, n, m*ra2(k));
            Bk = permute(Bk, [1,3,2,4]);
            Bk = reshape(Bk, rw1*n*rx1*m, ra2(k));
            tmp = reshape(WAX2{k}, ra2(k), rw2*rx2);
            Bk = Bk*tmp;
            Bk = reshape(Bk, rw1*n, rx1*m, rw2, rx2);
            Bk = permute(Bk, [1, 3, 2, 4]);
            Bk = reshape(Bk, rw1*n*rw2, rx1*m*rx2);
            B = B+Bk;
        end
    end
end
