% 例子1
filename = 'cdde3.rua';
A = hb_to_msm ( filename );
B = sum(A,2);
tol = 1e-6;
maxit = 10000;
%s = 1
[X1,r1] = IDR_(A,B,1,tol);
setup = struct('type','ilutp','droptol',1e-6);
semilogy(0:length(r1)-1,r1,'r-o');
hold on
%s = 2
[X2,r2] = IDR_(A,B,2,tol);
semilogy(0:length(r2)-1,r2,'b-o');
hold on
%s = 4
[X3,r3] = IDR_(A,B,4,tol);
semilogy(0:length(r3)-1,r3,'g-o');
hold on
%bicgstab
[X4,fl4,rr4,it4,rv4] = bicgstab(A,B,tol,maxit);
semilogy(0:length(rv4)-1,rv4/norm(B),'k-o')
hold on
yline(tol,'m--');
legend('s=1','s=2','s=4','bicgstab','tol');
xlabel('Iteration number')
ylabel('Relative residual')
%% 例子2
filename = 'cdde1.rua';
A = hb_to_msm ( filename );
B = sum(A,2);
tol = 1e-6;
maxit = 10000;
%s = 1
[X1,r1] = IDR_(A,B,1,tol);
figure
setup = struct('type','ilutp','droptol',1e-6);
semilogy(0:length(r1)-1,r1,'r-o');
hold on
%s = 2
[X2,r2] = IDR_(A,B,2,tol);
semilogy(0:length(r2)-1,r2,'b-o');
hold on
%s = 4
[X3,r3] = IDR_(A,B,4,tol);
semilogy(0:length(r3)-1,r3,'g-o');
hold on
%bicgstab
[X4,fl4,rr4,it4,rv4] = bicgstab(A,B,tol,maxit);
semilogy(0:length(rv4)-1,rv4/norm(B),'k-o')
hold on
yline(tol,'m--');
legend('s=1','s=2','s=4','bicgstab','tol');
xlabel('Iteration number')
ylabel('Relative residual')
%% 例子3
filename = 'UTM1700a.rua';
A = hb_to_msm ( filename );
B = sum(A,2);
tol = 1e-6;
maxit = 10000;
%s = 1
[X1,r1] = IDR_(A,B,1,tol);
figure
setup = struct('type','ilutp','droptol',1e-6);
semilogy(0:length(r1)-1,r1,'r-o');
hold on
%s = 2
[X2,r2] = IDR_(A,B,2,tol);
semilogy(0:length(r2)-1,r2,'b-o');
hold on
%s = 4
[X3,r3] = IDR_(A,B,4,tol);
semilogy(0:length(r3)-1,r3,'g-o');
hold on
%bicgstab
[X4,fl4,rr4,it4,rv4] = bicgstab(A,B,tol,maxit);
semilogy(0:length(rv4)-1,rv4/norm(B),'k-o')
hold on
yline(tol,'m--');
legend('s=1','s=2','s=4','bicgstab','tol');
xlabel('Iteration number')
ylabel('Relative residual')
%% 例子4
filename = 'fidap001.rua';
A = hb_to_msm ( filename );
B = sum(A,2);
tol = 1e-6;
maxit = 10000;
%s = 1
[X1,r1] = IDR_(A,B,1,tol);
figure
setup = struct('type','ilutp','droptol',1e-6);
semilogy(0:length(r1)-1,r1,'r-o');
hold on
%s = 2
[X2,r2] = IDR_(A,B,2,tol);
semilogy(0:length(r2)-1,r2,'b-o');
hold on
%s = 4
[X3,r3] = IDR_(A,B,4,tol);
semilogy(0:length(r3)-1,r3,'g-o');
hold on
%bicgstab
[X4,fl4,rr4,it4,rv4] = bicgstab(A,B,tol,maxit);
semilogy(0:length(rv4)-1,rv4/norm(B),'k-o')
hold on
yline(tol,'m--');
legend('s=1','s=2','s=4','bicgstab','tol');
xlabel('Iteration number')
ylabel('Relative residual')
%% 例子5
filename = 'fidap022.rua';
A = hb_to_msm ( filename );
B = sum(A,2);
tol = 1e-6;
maxit = 10000;
%s = 1
[X1,r1] = IDR_(A,B,1,tol);
figure
setup = struct('type','ilutp','droptol',1e-6);
semilogy(0:length(r1)-1,r1,'r-o');
hold on
%s = 2
[X2,r2] = IDR_(A,B,2,tol);
semilogy(0:length(r2)-1,r2,'b-o');
hold on
%s = 4
[X3,r3] = IDR_(A,B,4,tol);
semilogy(0:length(r3)-1,r3,'g-o');
hold on
%bicgstab
[X4,fl4,rr4,it4,rv4] = bicgstab(A,B,tol,maxit);
semilogy(0:length(rv4)-1,rv4/norm(B),'k-o')
hold on
yline(tol,'m--');
legend('s=1','s=2','s=4','bicgstab','tol');
xlabel('Iteration number')
ylabel('Relative residual')
%% 例子6
filename = 'sherman5.rua';
A = hb_to_msm ( filename );
B = sum(A,2);
tol = 1e-6;
maxit = 10000;
%s = 1
[X1,r1] = IDR_(A,B,1,tol);
figure
setup = struct('type','ilutp','droptol',1e-6);
semilogy(0:length(r1)-1,r1,'r-o');
hold on
%s = 2
[X2,r2] = IDR_(A,B,2,tol);
semilogy(0:length(r2)-1,r2,'b-o');
hold on
%s = 4
[X3,r3] = IDR_(A,B,4,tol);
semilogy(0:length(r3)-1,r3,'g-o');
hold on
%bicgstab
[X4,fl4,rr4,it4,rv4] = bicgstab(A,B,tol,maxit);
semilogy(0:length(rv4)-1,rv4/norm(B),'k-o')
hold on
yline(tol,'m--');
legend('s=1','s=2','s=4','bicgstab','tol');
xlabel('Iteration number')
ylabel('Relative residual')
%% 例子7
filename = 'saylr4.rua';
A = hb_to_msm ( filename );
B = sum(A,2);
tol = 1e-6;
maxit = 10000;
%s = 1
[X1,r1] = IDR_(A,B,1,tol);
figure
setup = struct('type','ilutp','droptol',1e-6);
semilogy(0:length(r1)-1,r1,'r-o');
hold on
%s = 2
[X2,r2] = IDR_(A,B,2,tol);
semilogy(0:length(r2)-1,r2,'b-o');
hold on
%s = 4
[X3,r3] = IDR_(A,B,4,tol);
semilogy(0:length(r3)-1,r3,'g-o');
hold on
%bicgstab
[X4,fl4,rr4,it4,rv4] = bicgstab(A,B,tol,maxit);
semilogy(0:length(rv4)-1,rv4/norm(B),'k-o')
hold on
yline(tol,'m--');
legend('s=1','s=2','s=4','bicgstab','tol');
xlabel('Iteration number')
ylabel('Relative residual')
%% 例子8
filename = 'sherman4.rua';
A = hb_to_msm ( filename );
B = sum(A,2);
tol = 1e-6;
maxit = 10000;
%s = 1
[X1,r1] = IDR_(A,B,1,tol);
figure
setup = struct('type','ilutp','droptol',1e-6);
semilogy(0:length(r1)-1,r1,'r-o');
hold on
%s = 2
[X2,r2] = IDR_(A,B,2,tol);
semilogy(0:length(r2)-1,r2,'b-o');
hold on
%s = 4
[X3,r3] = IDR_(A,B,4,tol);
semilogy(0:length(r3)-1,r3,'g-o');
hold on
%bicgstab
[X4,fl4,rr4,it4,rv4] = bicgstab(A,B,tol,maxit);
semilogy(0:length(rv4)-1,rv4/norm(B),'k-o')
hold on
yline(tol,'m--');
legend('s=1','s=2','s=4','bicgstab','tol');
xlabel('Iteration number')
ylabel('Relative residual')