function [x_ap, resvec_est, times] = sgmres(A, b, x0, d, l, B_mode, ...
    SE_mode, tol, QR_mode, s, spec_rec)
%sGMRES   sketched Generalized Minimum Residual Method.
%   x = sgmres(A, b, x0, d, l) attempts to solve the linear algebraic
%   system Ax = b, where A is nonsingular square N-by-N array, b is N-by-1 
%   array, x0 is N-by-1 array initial guess. Maximum number of iterations
%   d and truncation parameter l for incomplete Arnoldi process.
%
%   x = sgmres(A, b, x0, d, l, B_mode) specifies Krylov basis vector
%   generation. Default is "MGS" for (truncated) Arnoldi with MGS.
%   "Chebyshev" for Chebyshev polynomial basis (then l is unrelevant).
%
%   x = sgmres(A, b, x0, d, l, B_mode, SE_mode) specifies applyed  subspace 
%   embeddig. Default is "dct2" for subsampled random cosin transform. "fft"
%   uses subsampled random fourier transform. "Gauss" uses Gaussian
%   embedding and "id" no embedding at all.
%
%   x = sgmres(A, b, x0, d, l, B_mode, SE_mode, tol) specifies the
%   tolerance which must be reached to terminate before maximum iteration d.
%   Default is -inf.
%
%   x = sgmres(A, b, x0, d, l, B_mode, SE_mode, tol, QR_mode) specifies QR
%   factorization for solving for the final solution. Default is "thin" for
%   thin QR with MGS, which is faster but less stable. "full" uses
%   full QR with Givens rotations for updating, which is slower but more stable.
%
%   x = sgmres(A, b, x0, d, l, B_mode, SE_mode, tol, QR_mode, s) specifies
%   subspace embedding dimension. Default is s = 2*(d+1).
%   
%   x = sgmres(A, b, x0, d, l, B_mode, SE_mode, tol, QR_mode, s, spec_rec) 
%   necessary if B_mode = "Chebyshev". Specifies rectangle 
%   {z : c − a < Re(z) < c + a, c − b < Imag(z) < c + b} which encloses the 
%   spectrum of A. spec_rec is 1-by-3 array with spec_rec = [c, a, b].
%
%   [x, resvec_est] = sgmres(A, b, x0, d, l,...) also returns N-by-1 array
%   of estimates residual norms at each iteration.
%
%   [x, resvec_est, times] = sgmres(A, b, x0, d, l,...) also returns N-by-1 
%   array of time needed to every iteration.
%
%   Examples
%   x = sgmres(A, b, x0, d, 2, "MGS", "dct2", 1e-10, "full")
%   sGMRES with l = 2 truncated Arnoldi, SRCT subspace embedding and full
%   QR factorization with Givens rotations.
%
%   x = sgmres(A, b, x0, d, d, "MGS", "id") standard GMRES but slower due to
%   different QR factorization computation with maxiter d.

%   References
%   Yuji Nakatsukasa and Joel A. Tropp. “Fast & Accurate Randomized Algorithms
%   for Linear Systems and Eigenvalue Problems". In: arXiv e-prints (Oct. 2021). 
%   doi:10.48550/arXiv.2111.00113

tic;
[n,m] = size(A);
if (n ~= m) || (size(b,1) ~= n) || (size(x0,1) ~= n)
    error('Dimensions of linear system components does not match.')
end
if ~exist('s','var')
    % check if subspace embedding dimension is not specified
    s = 2*(d+1);
end

if ~exist('SE_mode','var')
    % check if Krylov basis mode is given and if not set it to MGS
    SE_mode = "dct2";
end

if ((SE_mode == "fft") || (SE_mode == "dct2")) && (s > floor(n/2))
    warning(['Subspace embedding dimension too great for meaningful sketching' ...
        'Setting max s = floor(n/2)'])
    s = floor(n/2);
end

if ~exist('tol','var')
    % check if tolerance is given and if not set it to -inf, to prevent 
    % stopping befor d itererations
    tol = -inf;
end

if ~exist('QR_mode','var')
    % check if QR mode is given and if not set it to thin QR-MGS
    QR_mode = "thin";
end

if ~exist('B_mode','var')
    % check if Krylov basis mode is given and if not set it to MGS
    B_mode = "MGS";
end

% initialize workspace
normb = norm(b);
times = zeros(d,1);
B = zeros(n,d); % Krylov Basis matrix
r = b-A*x0; % initial residual vector
[E, perm] = buildSubParam(n, s, SE_mode);

% initialize truncated Arnoldi/Chebyshev: build 1st basis vector
B(:,1) = r/norm(r);
g = applySubEmb(r, s, SE_mode, perm, E); % sketched initial residual
w = A*B(:,1);
if QR_mode == "thin"
    [Q,R] = qr(applySubEmb(w, s, SE_mode, perm, E),0);
elseif QR_mode == "full"
    [Q,R] = qr(applySubEmb(w, s, SE_mode, perm, E));
else
    error("Unkown QR update mode.")
end
r_est = g - (Q(:,1)'*g)*Q(:,1);
resvec_est = zeros(d,1);
resvec_est(1) = norm(r_est);
times(1) = toc;

% Build Krylov subspace Basis
if B_mode == "MGS"  % truncated Arnoldi with Modified Gram-Schmidt
    for j = 1:d-1
        for i = max(j-l+1,1):j
            w = w - B(:,i)*(B(:,i)'*w);
        end
        B(:,j+1) = w/norm(w);
        w = A*B(:,j+1);
        if QR_mode == "thin"
            [Q,R] = appendThinQR(Q,R,applySubEmb(w, s, SE_mode, perm, E));
        else
            [Q,R] = qrinsert(Q,R,j+1,applySubEmb(w, s, SE_mode, perm, E));
        end
        r_est = r_est - (Q(:,j+1)'*g)*Q(:,j+1);
        resvec_est(j+1) = norm(r_est);
        times(j+1) = toc;
        if resvec_est(j+1)/normb < tol % check for convergence
            B = B(:,1:j+1);
            resvec_est = resvec_est(1:j+1);
            times = times(1:j+1);
            break;
        end
    end

elseif B_mode == "Chebyshev"
    c = spec_rec(1);
    alpha = spec_rec(2);
    beta = spec_rec(3);
    rho = max(alpha,beta);   
    Ac = spdiags(spdiags(A,0)-c,0,A);
    B(:,2) = Ac*B(:,1)/(2*rho);
    B(:,3) = (Ac*B(:,2) - (alpha^2 - beta^2)*B(:,1)/(4*rho))/rho;
    for j = 2:d
        if j <= d-2
            B(:,j+2) = (Ac*B(:,j+1) - (alpha^2 - beta^2)*B(:,j)/(4*rho))/rho;
        end
        B(:,j) = B(:,j)/norm(B(:,j));
        if QR_mode == "thin"
            [Q,R] = appendThinQR(Q,R,applySubEmb(A*B(:,j), s, SE_mode, perm, E));
        else
            [Q,R] = qrinsert(Q,R,j,applySubEmb(A*B(:,j), s, SE_mode, perm, E));
        end
        r_est = r_est - (Q(:,j)'*g)*Q(:,j);
        resvec_est(j) = norm(r_est);
        times(j) = toc;
        if (resvec_est(j)/normb < tol) || isnan(resvec_est(j))
            B = B(:,1:j);
            resvec_est = resvec_est(1:j);
            times = times(1:j);
            break;
        end
    end
else
    error('Unkown Krylov subspace basis.')
end
% solve least-squares problem with QR factorization
if QR_mode == "full"
    R = R(1:size(R,2),:);
    Q = Q(:,1:size(R,2));
end
opts.UT = true;
y = linsolve(R,Q'*g,opts);
x_ap = x0 + B*y;
times(end) = toc;
end