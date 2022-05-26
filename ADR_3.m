% clc,clear
% % 
% % % n = 5;
% % % A = [ 1,   0,   0,   5,  0;
% % %       0,   2,   0,   0,  3;
% % %       1,   3,   4,   0,  0;
% % %       0,   0 ,  4,   3,  0;
% % %       0,   0 ,  4,   0,  1];
% % % A = sparse(A);
% rng('default')
% m = 2;
% s = 4;
% 
% [A, rows, cols, entries] = mmread('cdde1.mtx');
% n = rows;
% B = rand([n,m]);
% P = rand([n,s*m]);
% [I, R__] = ADR_6(A,m,rows,s,B,P);
function [I, R__] = ADR_3(A,m,n,s,B,P)

    Delta_X = zeros([n,s*m]);
    Delta_R = zeros([n,s*m]);
    X_ = zeros([n,m,18]);
    R__ = zeros([n,m,18]);

    % K = K1*K2, K2=I
    options.type = 'nofill';
    options.milu = 'off';
    [L,U] = ilu(A,options);
    K = L*U;
%     K_1 = K;

%     B = rand([n,m]);

    R__(:,:,1) = B - A*(L*U*X_(:,:,1));
%     P = rand([n,s*m]);

    % Initial R1,...,Rs+1(in papaer is R0,...Rs),
    % save in R_
    for I=1:s
        W = K\R__(:,:,I);
        V = A*W;
        w = trace((K\V)'/K*R__(:,:,I))/trace((K\V)'/K*V);
        Delta_X(:,(I-1)*m+1:I*m) = w*R__(:,:,I);
        Delta_R(:,(I-1)*m+1:I*m) = -w*V;
        X_(:,:,I+1) = X_(:,:,I) + Delta_X(:,(I-1)*m+1:I*m);
        R__(:,:,I+1) = R__(:,:,I) + Delta_R(:,(I-1)*m+1:I*m);
    end
    M = P'*Delta_R;
    h = P'*R__(:,:,s+1);
    j = 1;
    I = s+1;

    
    while biggest(R__(:,:,I),B,m) > 1e-8
        for k = 0:s
            C = M\h;
            Q = -Delta_R*C;
            V = R__(:,:,I) + Q;            
            if k == 0                    
                T = A/K*V;    
                w = trace((K\T)'/K*V)/trace((K\T)'/K*T);                   
                Delta_R(:,(j-1)*m+1:j*m) = Q - w*T;
                Delta_X(:,(j-1)*m+1:j*m) = -Delta_X*C + w*V;                
            else
                Delta_X(:,(j-1)*m+1:j*m) = -Delta_X*C + w*V;                  
                Delta_R(:,(j-1)*m+1:j*m) = -A/K*Delta_X(:,(j-1)*m+1:j*m);                
            end            
            R__(:,:,I+1) = R__(:,:,I) + Delta_R(:,(j-1)*m+1:j*m);
            X_(:,:,I+1) = X_(:,:,I) + Delta_X(:,(j-1)*m+1:j*m);
            Delta_m = P'*Delta_R(:,(j-1)*m+1:j*m);
            M(:,(j-1)*m+1:j*m) = Delta_m;
            h = h + Delta_m;
            I = I + 1;
            j = j + 1;
            j = mod((j-1),s) + 1;               
        end       
    end
    
%     X_(:,:,I) = K\X_(:,:,I);
end


function big = biggest(R,B,m)
    big = 0;
    for j = 1:m
        R_B = norm(R(:,j))/norm(B(:,j));
        if R_B > big
            big = R_B;
        end
    end
end






