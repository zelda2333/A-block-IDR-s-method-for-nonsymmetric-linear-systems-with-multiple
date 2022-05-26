 %function [X,r] = IDR_(A,B,s,tol)
%   filename = 'cdde3.rua';
%  filename = 'cdde1.rua';
%  filename = 'UTM1700a.rua';
%   filename = 'fidap001.rua';
%  filename = 'fidap022.rua';
%   filename = 'sherman5.rua';
%   filename = 'saylr4.rua';
%   filename = 'sherman4.rua';
% A = hb_to_msm ( filename );
% B = sum(A,2);
tol = 1e-6;
[A, rows, cols, entries] = mmread('cdde1.mtx');
n = rows;
%n = size(A,1);
% maxit = 10000;
m = 1;
% s = 1;
P = B;
% digits(3);
X = cell(1,100*n);
X(1) = {zeros(n,m)};
% vpa(cell2mat(X(1)));
R = cell(1,100*n);
R(1) = {B-A*cell2mat(X(1))};
% vpa(cell2mat(R(1)));
r = zeros(2*n+1,1);
r(1) = 1;
for i = 1:s
    V = A*cell2mat(R(i));
    omega = trace(V'*cell2mat(R(i)))/trace(V'*V);
    delta_X(:,(i-1)*m+1:i*m) = omega*cell2mat(R(i));
    delta_R(:,(i-1)*m+1:i*m) = -omega*V;
    X(i+1) = {cell2mat(X(i))+delta_X(:,(i-1)*m+1:i*m)};
%     vpa(cell2mat(X(i+1)));
    R(i+1) = {cell2mat(R(i))+delta_R(:,(i-1)*m+1:i*m)};
    r(i+1) = vecnorm(cell2mat(R(i+1)))./norm(B);
%     vpa(cell2mat(R(i+1)));
end
j = 2 ;
i = s + 1 ;
M = P'*delta_R;
h = P'*cell2mat(R(i));
% r = 1;
while r(i) > tol
    for k = 1:s+1
        C = M\h;
        Q = -delta_R*C;
        V = cell2mat(R(i))+Q;
        if k == 1
            T = A*V;
            omega = trace(T'*V)/trace(T'*T);
            delta_R(:,(j-2)*m+1:(j-1)*m) = Q-omega*T;
            delta_X(:,(j-2)*m+1:(j-1)*m) = -delta_X*C+omega*V;
        else
            delta_X(:,(j-2)*m+1:(j-1)*m) = -delta_X*C+omega*V;
            delta_R(:,(j-2)*m+1:(j-1)*m) = -A*delta_X(:,(j-2)*m+1:(j-1)*m);
        end
        R(i+1) = {cell2mat(R(i))+delta_R(:,(j-2)*m+1:(j-1)*m)};
        r(i+1) = vecnorm(cell2mat(R(i+1)))./vecnorm(B);
        %         vpa(cell2mat(R(i+1)));
        X(i+1) = {cell2mat(X(i))+delta_X(:,(j-2)*m+1:(j-1)*m)};
        %         vpa(cell2mat(X(i+1)));
        delta_m = P'*delta_R(:,(j-2)*m+1:(j-1)*m);
        M(:,(j-2)*m+1:(j-1)*m) = delta_m;
        h = h+delta_m;
        i = i+1;
        j = j+1;
        j = mod(j-2,s)+2;
        %         r(i+1) = vecnorm(cell2mat(R(i+1)))./vecnorm(B);
    end
end
