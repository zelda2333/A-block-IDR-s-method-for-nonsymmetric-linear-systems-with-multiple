clc,clear
rng('default')
m = 5;
s = 4;

[A, rows, cols, entries] = mmread('cdde1.mtx');
n = rows;
B = rand([n,m]);
P = rand([n,s*m]);
[i, R_] = ADR_5(A,m,rows,s,B,P);
function [i, R_] = ADR_5(A,m,n,s,B,P)
    Delta_X = zeros([n,s*m]);
    Delta_R = zeros([n,s*m]);
    X_ = zeros([n,m,18]);
    R_ = zeros([n,m,18]);

%     B = rand([n,m]);
    R_(:,:,1) = B - A*X_(:,:,1);
%     P = rand([n,s*m]);


    % Initial R1,...,Rs+1(in papaer is R0,...Rs),
    % save in R_
    for i=1:s
        V = A*R_(:,:,i);
        w = trace(V'*R_(:,:,i))/trace(V'*V);
        Delta_X(:,(i-1)*m+1:i*m) = w*R_(:,:,i);
        Delta_R(:,(i-1)*m+1:i*m) = -w*V;
        X_(:,:,i+1) = X_(:,:,i) + Delta_X(:,(i-1)*m+1:i*m);
        R_(:,:,i+1) = R_(:,:,i) + Delta_R(:,(i-1)*m+1:i*m);
    end
    
    j = 1;
    i = s+1;
    M = P'*Delta_R;
    h = P'*R_(:,:,i);


%     && i<2500
    while biggest(R_(:,:,i),B,m) > 1e-8
        for k = 0:s
            C = M\h;
            Q = -Delta_R*C;
            V = R_(:,:,i) + Q;
            if k == 0
                T = A*V;
                w = trace(T'*V)/trace(T'*T);
                Delta_R(:,(j-1)*m+1:j*m) = Q - w*T;
                Delta_X(:,(j-1)*m+1:j*m) = -Delta_X*C + w*V;
            else
                Delta_X(:,(j-1)*m+1:j*m) = -Delta_X*C + w*V;
                Delta_R(:,(j-1)*m+1:j*m) = -A*Delta_X(:,(j-1)*m+1:j*m);
            end
            R_(:,:,i+1) = R_(:,:,i) + Delta_R(:,(j-1)*m+1:j*m);
            X_(:,:,i+1) = X_(:,:,i) + Delta_X(:,(j-1)*m+1:j*m);
            Delta_m = P'*Delta_R(:,(j-1)*m+1:j*m);
            M(:,(j-1)*m+1:j*m) = Delta_m;
            h = h + Delta_m;
            i = i + 1;
            j = j + 1;
            j = mod((j-1),s) + 1;                 
        end
    end
end

function big = biggest(R,B,m)
    big = 0;
    for j = 1:m
        R_B = norm(R(:,j)) / norm(B(:,j));
        if R_B > big
            big = R_B;
        end
    end
end