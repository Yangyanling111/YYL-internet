function [X1,X2,X3] = Runge_Kutta1(h,N,f,g,z,tau,D,f0,g0,z0)
%%
%                                    %
% 计算带有时滞/随机/碰撞的微分方程组 %
%                                    %
%%
% h——步长
% N——计算次数
% f,g——待计算函数
% tau——时滞
% D——噪声强度
% f0、g0——初值
%%
GSN = [];
for i = 1:length(D)
    GsnTemp = wgn(1,N,D(i));
    GSN = [GSN;GsnTemp]; 
end
%%
X1 = ones(1,N)*f0;
X2 = ones(1,N)*g0;
X3 = ones(1,N)*z0;
%%
for i = 1:1:N-1
    if (i-1)*h <= tau
    %if i*h <= tau
        x1tau = 0;
    else
        x1tau = X1(ceil(i-(tau-mod(tau,h))/h));
    end
    x1 = X1(i);x2 = X2(i);x3 = X3(i);
    Kf1 = feval(f,x1tau,GSN(1,i),GSN(2,i),x1,x2,x3);
    Kg1 = feval(g,x1tau,GSN(1,i),GSN(2,i),x1,x2,x3);
    Kz1 = feval(z,x1tau,GSN(1,i),GSN(2,i),x1,x2,x3);
    Kf2 = feval(f,x1tau,GSN(1,i),GSN(2,i),x1+Kf1*h/2,x2+Kg1*h/2,x3+Kz1*h/2);
    Kg2 = feval(g,x1tau,GSN(1,i),GSN(2,i),x1+Kf1*h/2,x2+Kg1*h/2,x3+Kz1*h/2);
    Kz2 = feval(z,x1tau,GSN(1,i),GSN(2,i),x1+Kf1*h/2,x2+Kg1*h/2,x3+Kz1*h/2);
    Kf3 = feval(f,x1tau,GSN(1,i),GSN(2,i),x1+Kf2*h/2,x2+Kg2*h/2,x3+Kz2*h/2);
    Kg3 = feval(g,x1tau,GSN(1,i),GSN(2,i),x1+Kf2*h/2,x2+Kg2*h/2,x3+Kz2*h/2);
    Kz3 = feval(z,x1tau,GSN(1,i),GSN(2,i),x1+Kf2*h/2,x2+Kg2*h/2,x3+Kz2*h/2);
    Kf4 = feval(f,x1tau,GSN(1,i),GSN(2,i),x1+Kf3*h,x2+Kg3*h,x3+Kz3*h);
    Kg4 = feval(g,x1tau,GSN(1,i),GSN(2,i),x1+Kf3*h,x2+Kg3*h,x3+Kz3*h);
    Kz4 = feval(z,x1tau,GSN(1,i),GSN(2,i),x1+Kf3*h,x2+Kg3*h,x3+Kz3*h);

    X1(i+1) = X1(i) + (Kf1 + 2*Kf2 + 2*Kf3 + Kf4)*h / 6;
    X2(i+1) = X2(i) + (Kg1 + 2*Kg2 + 2*Kg3 + Kg4)*h / 6;
    X3(i+1)=X3(i) + (Kz1 + 2*Kz2 + 2*Kz3 + Kz4)*h / 6;
end