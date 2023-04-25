clc
clear
clear all
global mu1 mu2 mu3 G
mu1=0.0545413;%0.01  0.0324413 0.0545413
mu2=-0.005;
mu3=0.0325413;
zz=10;
z=-zz:0.01:zz;
[z1,z2]=meshgrid(z);
A=sqrt(z1.^2+z2.^2);
C=2;
p=C./mu3.*A.^(2.*mu1./mu3-2).*exp(mu2.*A.^2./mu3);
pst=real(p);
meshz(z1,z2,pst)
h1=xlabel('$z_1$','FontSize',20);
 h2=ylabel('$z_2$','FontSize',20);
h3=zlabel('$p_{st}$','FontSize',20);
set(gca,'FontSize',20,'Fontname', 'Times New Roman');
 set(h1,'Interpreter','latex');
 set(h2,'Interpreter','latex');
 set(h3,'Interpreter','latex');
mesh()
