% ---------------------------------------------------------------------------------
% FUNCTION INFORMATION (c) 2023 Telecommunications Circuits Laboratory, EPFL
% ---------------------------------------------------------------------------------
% name  : BinaryGaussianEliminationForSCLDPC
% descr : perform Gaussian Elimination to transform H as a standard format [P|I]

function [H, Identity, IdentityInv] = BinaryGaussianEliminationForSCLDPC(H, M, N, w)
K           = N-M;
Identity    = eye(M);
IdentityInv = eye(M);

% move last Mb to the beginning
H(:, :, 1) = [H(:, K+1:end, 1) H(:, 1:K, 1)];

if isequal(H(:, 1:M, 1), eye(M))
    fprintf('Orignal matrix has already satisfied a standard format.\n');
else
    fprintf('Perform Gaussian elimination by row transformation to get a standard format.\n');

    % start to do Gaussian Elimination
    for i = 1:M
        pivot = -1;
        for j = i:M
            if H(j, i, 1) == 1
                pivot = j;
                break;
            end
        end

        if (pivot ~= -1)
            % Add
            if i ~= pivot
                for k = 1:N
                    for i_w = 1:w
                        H(i, k, i_w) = mod(H(i, k, i_w) + H(pivot, k, i_w), 2);
                    end
                end

                for k = 1:M
                    Identity(i, k)        = mod(Identity(i, k) + Identity(pivot, k), 2);
                    IdentityInv(k, pivot) = mod(IdentityInv(k, pivot) + IdentityInv(k, i), 2);
                end
            end

            % Eliminate other rows (from top to bottom)
            for j = i+1:M
                if H(j, i, 1) == 1
                    for k = 1:N
                        for i_w = 1:w
                            H(j, k, i_w) = mod(H(j, k, i_w) + H(i, k, i_w), 2);
                        end
                    end

                    for k = 1:M
                        Identity(j, k)    = mod(Identity(j, k) + Identity(i, k), 2);
                        IdentityInv(k, i) = mod(IdentityInv(k, i) + IdentityInv(k, j), 2);
                    end
                end
            end
        else
            fprintf("Can not find this row %0d\n", i);
        end
    end

    % Eliminate other rows (from bottom to top)
    for i = M:-1:1
        for j = i-1:-1:1
            if H(j, i) == 1
                for k = 1:N
                    for i_w = 1:w
                        H(j, k, i_w) = mod(H(j, k, i_w) + H(i, k, i_w), 2);
                    end
                end

                for k = 1:M
                    Identity(j, k)    = mod(Identity(j, k) + Identity(i, k), 2);
                    IdentityInv(k, i) = mod(IdentityInv(k, i) + IdentityInv(k, j), 2);
                end
            end
        end
    end

end

% move back
H(:, :, 1) = [H(:, M+1:end, 1) H(:, 1:M, 1)];
end
