function [E, perm] = buildSubParam(n, s, SE_mode)
% Builds parameters for subspace embedding

if SE_mode == "fft"
    phi = rand(n,1)*2*pi; % Uniform random variable sample form (0,2*pi).
    E = cos(phi) + 1i*sin(phi); % Steinhaus random variable samples.
    perm = datasample(1:n,s);
elseif SE_mode == "dct2"
    E = binornd(1,0.5,n,1)*2-1; % Rademacher random variable samples.
    perm = datasample(1:n,s);
elseif SE_mode == "id" % Dummy parameters for no dimension reduction.
    E = 1;
    perm = 1;
elseif SE_mode == "Gauss"
    E = normrnd(0,1,s,n); % Standard normal random variable samples.
    perm = 1;
else
    warning('Unkown subspace embedding method.')
end
end