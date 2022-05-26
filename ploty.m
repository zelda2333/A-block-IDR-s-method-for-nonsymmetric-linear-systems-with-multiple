clear,clc
LW = 'LineWidth'; lw = 1;
% FS = 'FontSize'; fs = 20;
MS = 'MarkerSize'; ms = 12;
IN = 'interpret';
LX = 'latex';

rng('default')
m = 20;
s = 4;

[A, rows, cols, entries] = mmread('cdde1.mtx');
n = rows;
B = rand([n,m]);
P = rand([n,s*m]);


% ADR_1
iter_ = 0;
r2 = zeros([n,1]);
norm_r2 = zeros([m*s,1]);

tic
for i = 1:m
    b = B(:,i);
    p = P(:,s*(i-1)+1 : s*i);
    [i_, r_] = ADR_1(A,n,s,b,p);
    iter_ = iter_ + i_;
    r2(:,(iter_-i_)+1 : iter_) = r_;    
end
toc


for j = 1:iter_
    norm_r2(j) = norm(r2(:,j));
end

r2__log = log10(norm_r2);

% ADR_2  
% [i, R_] = ADR_2(A,m,n,s,B,P);
% max_r2 = zeros([i,1]);
% for iter = 1:i
%    max_r2(iter) = biggest(R_(:,:,iter),m);
% end
% r2_log = log10(max_r2);


% ADR_3 
tic
[I, R__] = ADR_3(A,m,n,s,B,P);
toc

max_R2 = zeros([I,1]);
norm_max_R2 = zeros([I,1]);
for iter = 1:I
   max_R2(iter) = biggest(R__(:,:,iter),m);
   norm_max_R2(iter) = norm(max_R2(iter));
end
R2_log = log10(norm_max_R2);

% % total
% figure(1)
% clf
% x = 5 : 5: 40;
% IDR = [270, 530, 800, 1095, 1335, 1615, 1915, 2140];
% block_IDR = [150, 250, 375,500, 625, 600, 875, 1000];
% % block_bigcstab = [0.4, 0.45, 0.6, 0.9, 1, 1.05, 1.25, 1.3];
% % plot(x, IDR, 'bs--',x, block_IDR, 'r*--',x,block_bigcstab,'ko-')
% plot(x, IDR, 'bs--',x, block_IDR, 'r*--')
% xlim([5 40])
% ylim([0 2500])
% xlabel('Numner of the right-hand sides')
% ylabel('Numner of the matrix-vector products')
% % legend('IDR(4)','block IDR(4)','block Bicgstab')
% legend('IDR(4)','block IDR(4)')
% legend('Location','northwest')
% set(gca,'FontSize',18)
% set(gca,'FontName','times')
% set(gca,'YTick',0:500:2500);
% set(gcf,'color',[1,1,1]);
% 
% % average 
% figure(2)
% clf
% x = 5 : 5: 40;
% IDR = [54, 53, 40, 54.75, 53.4, 53.83, 54.71, 53.5];
% block_IDR = [30, 25, 25,25, 25, 20, 25, 25];
% % block_bigcstab = [0.4, 0.45, 0.6, 0.9, 1, 1.05, 1.25, 1.3];
% % plot(x, IDR, 'bs--',x, block_IDR, 'r*--',x,block_bigcstab,'ko-')
% plot(x, IDR, 'bs--',x, block_IDR, 'r*--')
% xlim([5 40])
% ylim([0 60])
% xlabel('Numner of the right-hand sides')
% ylabel('Numner of the matrix-vector products')
% % legend('IDR(4)','block IDR(4)','block Bicgstab')
% legend('IDR(4)','block IDR(4)')
% legend('Location','southeast')
% set(gca,'FontSize',18)
% set(gca,'FontName','times')
% set(gca,'YTick',0:10:60);
% set(gcf,'color',[1,1,1]);
% 
% 
% % total
% figure(3)
% clf
% x = 5 : 5: 40;
% IDR = [13.78, 28.32, 41.37, 57.12, 68.88, 84.03, 99.33, 112.6];
% IDR_2_ = [0.03, 0.05, 0.06, 1.47, 0.79, 0.42, 1.07, 1.3];
% block_IDR = [1.50, 1.37, 1.32,1.46, 1.41, 1.13, 1.63, 1.56];
% % block_bigcstab = [0.4, 0.45, 0.6, 0.9, 1, 1.05, 1.25, 1.3];
% % plot(x, IDR, 'bs--',x, block_IDR, 'r*--',x,block_bigcstab,'ko-')
% plot(x, IDR, 'bs--',x, block_IDR, 'r*--')
% xlim([5 40])
% ylim([0 120])
% xlabel('Numner of the right-hand sides')
% ylabel('CPU time(second)')
% % legend('IDR(4)','block IDR(4)','block Bicgstab')
% legend('IDR(4)','block IDR(4)')
% legend('Location','northwest')
% set(gca,'FontSize',18)
% set(gca,'FontName','times')
% set(gca,'YTick',0:20:120);
% set(gcf,'color',[1,1,1]);
% 
% % average 
% figure(4)
% clf
% x = 5 : 5: 40;
% IDR = [13.78/5, 28.32/10, 41.37/15, 57.12/20, 68.88/25, 84.03/30, 99.33/36, 112.6/40];
% % IDR_2_ = [0.03, 0.05, 0.06, 1.47, 0.79, 0.42, 1.07, 1.3];
% block_IDR = [1.5/5, 1.37/10, 1.32/15,1.46/20, 1.41/25, 1.13/30, 1.63/36, 1.56/40];
% % block_bigcstab = [0.4, 0.45, 0.6, 0.9, 1, 1.05, 1.25, 1.3];
% % plot(x, IDR, 'bs--',x, block_IDR, 'r*--',x,block_bigcstab,'ko-')
% plot(x, IDR, 'bs--',x, block_IDR, 'r*--')
% xlim([5 40])
% ylim([0 3])
% xlabel('Numner of the right-hand sides')
% ylabel('CPU time(second)')
% % legend('IDR(4)','block IDR(4)','block Bicgstab')
% legend('IDR(4)','block IDR(4)')
% legend('Location','southeast')
% set(gca,'FontSize',18)
% set(gca,'FontName','times')
% set(gca,'YTick',0:0.5:3);
% set(gcf,'color',[1,1,1]);

% m=10 
figure(5)
clf
x = 1:1:iter_;
X = 1:10:I*10;
IDR = r2__log;
block_IDR = R2_log;
%block_bigcstab = [0.4, 0.45, 0.6, 0.9, 1, 1.05, 1.25, 1.3];
%plot(x, IDR, 'bs--',x, block_IDR, 'r*--',x,block_bigcstab,'ko-')
plot(X, block_IDR, 'r*--',x, IDR, 'b--')
xlim([0 600])
ylim([-10 5])
xlabel('Numner of matrix-vector products')
ylabel('Maximum relative residual 2-nrom')
% legend('block IDR(4)','block Bicgstab','IDR(4)')
legend('block IDR(4)','IDR(4)')
legend('Location','southeast')
set(gca,'FontSize',18)
set(gca,'FontName','times')
set(gca,'YTick',-10:5:5);
set(gca,'XTick',0:100:600);
set(gcf,'color',[1,1,1]);

% m=20 
figure(5)
clf
x = 1:1:iter_;
X = 1:20:I*20;
IDR = r2__log;
block_IDR = R2_log;
%block_bigcstab = [0.4, 0.45, 0.6, 0.9, 1, 1.05, 1.25, 1.3];
%plot(x, IDR, 'bs--',x, block_IDR, 'r*--',x,block_bigcstab,'ko-')
plot(X, block_IDR, 'r*--',x, IDR, 'b--')
xlim([0 1200])
ylim([-10 5])
xlabel('Numner of matrix-vector products')
ylabel('Maximum relative residual 2-nrom')
% legend('block IDR(4)','block Bicgstab','IDR(4)')
legend('block IDR(4)','IDR(4)')
legend('Location','southeast')
set(gca,'FontSize',18)
set(gca,'FontName','times')
set(gca,'YTick',-10:5:5);
set(gca,'XTick',0:200:1200);
set(gcf,'color',[1,1,1]);

function big = biggest(R,m)
    big = 0;
    for j = 1:m
        R_norm = norm(R(:,j));
        if R_norm > big
            big = R_norm;
        end
    end
end