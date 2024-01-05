function C = applySubEmb(M, s, SE_mode, perm, E)
% applies implicitly a subspace embedding S = sqrt(n/s)*D*F*E to an 1/2D 
% array M (C = S*M). Assumes diagonal matrix is given as vector E with 
% Steinhaus/ Rademacher RV entries. Applies fft/dct2 with MATLAB build-in
% functions. Applies row sampling matrix D through row selection with 
% permutation array 1xn perm.

n = size(M,1);
if SE_mode == "fft"
    C = fft(E.*M);
    C = C(perm,:)/sqrt(s);
elseif SE_mode == "dct2"
    C = dct(E.*M,[],1,'Type',2);
    C = sqrt(n/s)*C(perm,:);
elseif SE_mode == "id"
    C = M;
elseif SE_mode == "Gauss"
    C = E*M/sqrt(s);
else
    warning('Unkown subspace embedding method.')
end
end