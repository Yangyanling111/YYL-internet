function run()
global beta
% clear all;
% clc
% close all
%%
k=0.0002;epsilon = 0.000001;c = 50;T=0.6;gam=0.3;omegaheng=0.001;tau=0.5;pxing = 1/c;beta=1.9834;
%tau=0.5,beta在1.94-1.9822之间
%beta=0.8,tau在1.2-1.25之间
%%
f0 = 0.1;
g0 = 0.1;
z0=0.1;
h = 0.01;
N = 40001;
K=0.282;
D = [2*K*pi,1];
%t = 0:h:(N-1)*h;
%%
% f -- x1'
% g -- x2'
f = @(x1tau,GSN1,GSN2,x1,x2,x3)(k*x1*(1/x2-c)+beta*(x1-x1tau)+epsilon*gam*(x1-pxing)*cos(omegaheng*x3)+epsilon^0.5*(x1-pxing)*GSN1);
g = @(x1tau,GSN1,GSN2,x1,x2,x3)((x1-x2)/T);
z = @(x1tau,GSN1,GSN2,x1,x2,x3)(1);
[X1,X2,X3] = Runge_Kutta(h,N,f,g,z,tau,D,f0,g0,z0); 
%%
t = 0:h:(N-1)*h;
% Drawing
%Time History
figure(1);
plot(t,X1);
hold on
h1=xlabel('$t$','FontSize',20);
h2=ylabel('$p$','FontSize',20);
set(gca,'FontSize',20,'Fontname', 'Times New Roman');
set(h1,'Interpreter','latex');
set(h2,'Interpreter','latex');

%Phase
figure(2)
plot(X1,X2);
h3=xlabel('$p$','FontSize',20);
h4=ylabel('$y$','FontSize',20);
set(gca,'FontSize',20,'Fontname', 'Times New Roman');
set(h3,'Interpreter','latex');
set(h4,'Interpreter','latex');

%probability density curves
figure(3)
ksdensity(X1);
h5=xlabel('$p$','FontSize',20);
h6=ylabel('$Density$','FontSize',20);
set(gca,'FontSize',20,'Fontname', 'Times New Roman');
set(h5,'Interpreter','latex');
set(h6,'Interpreter','latex');
