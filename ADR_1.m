% clc;
% clear all;
% function x = IDR(s,A)
% end
% s = 4;
% n = 5;
% A=[1,2,3,4,5;
%    0,1,2,3,4;
%    2,5,6,7,0;
%    4,3,0,5,2;
%    0,0,4,7,3];
% A = [ 1,   0,   0,   5,  0;
%       0,   2,   0,   0,  3;
%       1,   3,   4,   0,  0;
%       0,   0 ,  4,   3,  0;
%       0,   0 ,  4,   0,  1];
% [A, rows, cols, entries] = mmread('cdde1.mtx');
% n = rows;
% [i_, r_] = ADR_4(A,n,s);
% function [i_, r_] = ADR_4(A,n,s)
function [i_, r_] = ADR_1(A,n,s,b,P)
    Delta_X = zeros([n,s]);
    Delta_R = zeros([n,s]);
    x_ = zeros([n,3]);
    r_ = zeros([n,3]);
    options.type = 'nofill';
    options.milu = 'off';
    [L,U] = ilu(A,options);
    K = L*U;
%     b = rand([n,1]);
    % b = [0.908745548779810;0.785250310098426;0.628586468419828;0.643045675988833;0.205787129849464];
    r_(:,1) = b - A*K*x_(:,1);
%     P = rand([n,s]);
%     P = [0.628489999070803,0.227279660820927,0.853087124286285,0.521151453400346;0.472595799615831,0.179744729891939,0.739889715258732,0.102322127726926;0.539610163192950,0.432779310960602,0.551257461382844,0.249351836527347;0.251666841925854,0.465537685184856,0.682490882499042,0.245092448337264;0.0273053923262909,0.576615246639198,0.911462702037672,0.403972205512039];

    % Initial r1,...,rs+1(in papaer is r0,...rs),
    % save in r_
    for i_=1:s
        W = K\r_(:,i_);
        v = A*W;
%         v = A*r_(:,i_);
        w = trace((K\v)'/K*r_(:,i_))/trace((K\v)'/K*v);
        Delta_X(:,i_) = w*r_(:,i_);
        Delta_R(:,i_) = -w*v;
        x_(:,i_+1) = x_(:,i_) + Delta_X(:,i_);
        r_(:,i_+1) = r_(:,i_) + Delta_R(:,i_);   
    end
    M = P'*Delta_R;
    h = P'*r_(:,s+1);
    j = 1;
    i_ = s+1;
    % while norm(r_(:,i))/norm(b) > 1e-8

    while biggest(r_,b,i_)> 1e-8
        for k = 0:s
            c = M\h;
            q = -Delta_R*c;
            v = r_(:,i_) + q;
            % (Enter the next subspace Gl)
            if k == 0
                t = A/K*v;
                w = trace((K\t)'/K*v)/trace((K\t)'/K*t);  
                Delta_R(:,j) = q - w*t;
                Delta_X(:,j) = -Delta_X*c + w*v;
            else
                Delta_X(:,j) = -Delta_X*c + w*v;
                Delta_R(:,j) = -A/K*Delta_X(:,j);
            end
            % (Update approximate solutions xi)
            r_(:,i_+1) = r_(:,i_) + Delta_R(:,j);
            x_(:,i_+1) = x_(:,i_) + Delta_X(:,j);
            delta_m = P'*Delta_R(:,j);
            M(:,j) = delta_m;
            h = h + delta_m;
            i_ = i_ + 1;
            j = j + 1;
            j = mod((j-1),s) + 1;
        end
    end
end




function big = biggest(r_,b,i)
    big = 0;    
    R_B =norm(r_(:,i))/norm(b);
    if R_B > big
        big = R_B;
    end
    
end

