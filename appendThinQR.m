function [Q,R] = appendThinQR(Q,R,v)
% Updates a thin QR factoization A=Q1*R1 to [A,v]= Q2*R2 after the column 
% v is inserted as last column into A. Uses modified Gram-Schmidt.

%init
m = size(Q,2);
R=[R zeros(m,1)];
R=[R;zeros(1,m+1)];

for j = 1:m
    R(j,m+1) = Q(:,j).'*v;
    v = v - R(j,m+1)*Q(:,j);
end
R(m+1,m+1) = norm(v);
Q=[Q v/R(m+1,m+1)];
end