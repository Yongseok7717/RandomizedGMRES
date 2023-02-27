load('LS89.mat');
n=length(A);
b=load('RHS.dat');
tol=1e-10;
t=1000;
theta=Rademacher_r(t,n);
m=60;

[V,H,Z]=inner_outer_Arnoldi_RGS_ras(A,b,60,30,theta,6,6);
[Vc,Hc,Zc]=CGS_ras(A,b,60,30,6,6);
[Vm,Hm,Zm]=MGS_ras(A,b,60,30,6,6);

cond_r=zeros(m+1,1);
cond_m=zeros(m+1,1);
cond_c=zeros(m+1,1);

error_r=zeros(m,1);
error_m=zeros(m,1);
error_c=zeros(m,1);

cond_r(1)=cond(V(:,1));
cond_m(1)=cond(Vm(:,1));
cond_c(1)=cond(Vc(:,1));
    
for i=1:m
    cond_r(i+1)=cond(V(:,1:i+1));
    error_r(i)=norm(A*Z(:,1:i)-V(:,1:i+1)*H(1:i+1,1:i))/norm(A*Z(:,1:i));


    cond_m(i+1)=cond(Vm(:,1:i+1));
    error_m(i)=norm(A*Zm(:,1:i)-Vm(:,1:i+1)*Hm(1:i+1,1:i))/norm(A*Zm(:,1:i));


    cond_c(i+1)=cond(Vc(:,1:i+1));
    error_c(i)=norm(A*Zc(:,1:i)-Vc(:,1:i+1)*Hc(1:i+1,1:i))/norm(A*Zc(:,1:i));
end

% Stability figure
figure(1)
iter=1:m+1;
semilogy(iter,cond_c,'-*',iter,cond_m,'-o',...
    iter,cond_r,'-s');
xlim([0,m]);
legend('CGS','MGS','RGS, t=1000', 'Location','best')
title('Condition numbers of $$V_i$$','interpreter','latex' )
xlabel('iterations $$i$$','interpreter','latex')

% Approximation errors figure
figure(2)
iter=1:m;
semilogy(iter,error_c,'-*',iter,error_m,'-o',iter,error_r,'-s');
xlim([0,m])
legend('CGS','MGS','RGS, t=1000', 'Location','NorthWest')
title('Error $$|| AZ_i-V_{i+1}H_i||/|| AZ_i||$$','interpreter','latex')
xlabel('iterations $$i$$','interpreter','latex')

%% Compare CGS- and RGS-FGMRES-(M)DR
% CGS
[x_r,res_r,iter_r]=FGMRES_dr_ras(A,b,tol,60,30,0,6,12);
[x_dr,res_dr,iter_dr]=FGMRES_dr_ras(A,b,tol,60,30,30,6,12);
[x_mdr,res_mdr,iter_mdr]=FGMRES_mdr_ras(A,b,tol,60,30,30,6,12);

% RGS
[x_Rr,res_Rr,iter_Rr]=RGS_FGMRES_dr_ras(A,b,tol,60,30,0,1000,6,12);
[x_Rdr,res_Rdr,iter_Rdr]=RGS_FGMRES_dr_ras(A,b,tol,60,30,30,1000,6,12);
[x_Rmdr,res_Rmdr,iter_Rmdr]=RGS_FGMRES_mdr_ras(A,b,tol,60,30,30,1000,6,12);

figure(3)
semilogy(iter_r,res_r,'-',iter_dr,res_dr,':',...
    iter_mdr,res_mdr,'-.',iter_Rr,res_Rr,'-x',...
    iter_Rdr,res_Rdr,'-s',iter_Rmdr,res_Rmdr,'-o')
% 
legend('FGMRES(60,30)','FGMRES DR(60,30,30)','FGMRES MDR(60,30,30)',...
'RGS FGMRES(60,30)','RGS FGMRES DR(60,30,30)','RGS FGMRES MDR(60,30,30)')
xlim([0,120])
xlabel('Iterations')
ylabel('Relative residual norms')