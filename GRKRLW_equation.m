clc;clear;
close all;
% fourth-order compact scheme for generalized Rosenau-Kawahara-RLW equation
% Example 2
T=10; hh=0.1;tau=hh^2;
xl=-40; xr=60;  
N=round(T/tau); J=(xr-xl)/hh;
aa=1; bb=1;cc=2;dd=1;ee=1;ff=1;pp=4;
for j=1:J+1
    x(j)=xl+(j-1)*hh;
end
for kk=1:N+1
    t(kk)=(kk-1)*tau;
end
% exact solution ====================================================
st1=sqrt(127)-10;
st2=10*sqrt(127)-109;
st3=118-10*sqrt(127);
for j=1:J+1
    for kk=1:N+1
        uu(j,kk)=(40*(st1)^2/(3*st2))^(1/4)*sech((sqrt(st1))/3*...
            (x(j)-st3*t(kk)/st2));
    end
end
% numerical solution
% initial condition
for j=1:J+1
    u(j,1)=uu(j,1);
end
% boundary condition
for kk=2:N+1
    u(1,kk)=0;    u(2,kk)=0;    u(3,kk)=0;
    u(J-1,kk)=0;  u(J,kk)=0;    u(J+1,kk)=0;
end

% numerical solution====================go===================================
tic;
% step 1: u(:,2)==================================================
s1=cc+aa*hh^2/6;       s2=dd-hh^2/3;
s3=ee-dd*hh^2/4;       s4=ff-cc*hh^2/12;
A=zeros(J-5,J-5);B=A;

ss1=aa*tau/(4*hh);         ss2=s1*tau/(4*hh^3);
ss3=s2/hh^2;               ss4=s3/hh^4;
ss5=s4*tau/(4*hh^5);       ss6=ee/(6*hh^4);
ss7=bb*tau/(2*(pp+1)*hh);  ss8=bb*tau/(12*(pp+1)*hh);

hy1=ss5+ss6;
hy2=-ss2+ss4-4*ss5-6*ss6;
hy3=-ss1+2*ss2-ss3-4*ss4+5*ss5+15*ss6;
hy4=1+2*ss3+6*ss4-20*ss6;
hy5=ss1-2*ss2-ss3-4*ss4-5*ss5+15*ss6;
hy6=ss2+ss4+4*ss5-6*ss6;
hy7=-ss5+ss6;
for j=1:J-5
    hy8(j)=hy6+ss8*(u(j+3-2,1))^pp;
    hy9(j)=hy5+(ss7-2*ss8)*(u(j+3-1,1))^pp;
    hy10(j)=hy3-(ss7-2*ss8)*(u(j+3+1,1))^pp;
    hy11(j)=hy2-ss8*(u(j+3+2,1))^pp;
end

for j=1:J-5
    A(j,j)=hy4;     B(j,j)=hy4;
end
for j=1:J-6
    A(j,j+1)=hy5;   B(j,j+1)=hy10(j);
end
for j=1:J-7
    A(j,j+2)=hy6;   B(j,j+2)=hy11(j);
end
for j=1:J-8
    A(j,j+3)=hy7;   B(j,j+3)=hy1;
end
for j=2:J-5
    A(j,j-1)=hy3;   B(j,j-1)=hy9(j);
end
for j=3:J-5
    A(j,j-2)=hy2;   B(j,j-2)=hy8(j);
end
for j=4:J-5
    A(j,j-3)=hy1;   B(j,j-3)=hy7;
end
for i=1:J-5
    HH(i)=0;
    for j=1:J-5
         HH(i)=HH(i)+B(i,j)*u(j+3,1);
    end
    h(i)=HH(i);
end
%  solve a sevendiaginal system Ax=f  ==============================
c=diag(A);
d=[diag(A,1);0];
e=[diag(A,2);0;0];
f=[diag(A,3);0;0;0];
b=[0;diag(A,-1)];
a=[0;0;diag(A,-2)];
g=[0;0;0;diag(A,-3)];

n=length(h);

alpha(1)=c(1);          
beta(1)=d(1)/alpha(1);
q(1)=e(1)/alpha(1); 
k(1)=f(1)/alpha(1);

gama(2)=b(2);
alpha(2)=c(2)-gama(2)*beta(1);
beta(2)=(d(2)-gama(2)*q(1))/alpha(2);
q(2)=(e(2)-gama(2)*k(1))/alpha(2);
k(2)=f(2)/alpha(2);

z(3)=a(3);
gama(3)=b(3)-z(3)*beta(1);
alpha(3)=c(3)-z(3)*q(1)-gama(3)*beta(2);
beta(3)=(d(3)-z(3)*k(1)-gama(3)*q(2))/alpha(3);
q(3)=(e(3)-gama(3)*k(2))/alpha(3);
k(3)=f(3)/alpha(3);

for i=4:n
    m(i)=g(i);
end

for i=4:n-3
    z(i)=a(i)-m(i)*beta(i-3);
    gama(i)=b(i)-m(i)*q(i-3)-z(i)*beta(i-2);
    alpha(i)=c(i)-m(i)*k(i-3)-z(i)*q(i-2)-gama(i)*beta(i-1);
    beta(i)=(d(i)-z(i)*k(i-2)-gama(i)*q(i-1))/alpha(i);
    q(i)=(e(i)-gama(i)*k(i-1))/alpha(i);
    k(i)=f(i)/alpha(i);
end

z(n-2)=a(n-2)-m(n-2)*beta(n-5);
gama(n-2)=b(n-2)-m(n-2)*q(n-5)-z(n-2)*beta(n-4);
alpha(n-2)=c(n-2)-m(n-2)*k(n-5)-z(n-2)*q(n-4)-gama(n-2)*beta(n-3);
beta(n-2)=(d(n-2)-z(n-2)*k(n-4)-gama(n-2)*q(n-3))/alpha(n-2);
q(n-2)=(e(n-2)-gama(n-2)*k(n-3))/alpha(n-2);

z(n-1)=a(n-1)-m(n-1)*beta(n-4);
gama(n-1)=b(n-1)-m(n-1)*q(n-4)-z(n-1)*beta(n-3);
alpha(n-1)=c(n-1)-m(n-1)*k(n-4)-z(n-1)*q(n-3)-gama(n-1)*beta(n-2);
beta(n-1)=(d(n-1)-z(n-1)*k(n-3)-gama(n-1)*q(n-2))/alpha(n-1);

z(n)=a(n)-m(n)*beta(n-3);
gama(n)=b(n)-m(n)*q(n-3)-z(n)*beta(n-2);
alpha(n)=c(n)-m(n)*k(n-3)-z(n)*q(n-2)-gama(n)*beta(n-1);

% Ly=H
y(1)=h(1)/alpha(1);
y(2)=(h(2)-gama(2)*y(1))/alpha(2);
y(3)=(h(3)-z(3)*y(1)-gama(3)*y(2))/alpha(3);  
for i=4:n
    y(i)=(h(i)-m(i)*y(i-3)-z(i)*y(i-2)-gama(i)*y(i-1))/alpha(i);
end

% Ux=y
as(n)=y(n);
as(n-1)=y(n-1)-beta(n-1)*as(n);
as(n-2)=y(n-2)-beta(n-2)*as(n-1)-q(n-2)*as(n);
for i=n-3:-1:1
    as(i)=y(i)-beta(i)*as(i+1)-q(i)*as(i+2)-k(i)*as(i+3);
end
% end sovle seven diagnoal ===================================
for j=4:J-2
    u(j,2)=as(j-3);
end
% end step 1 ==========================================================

%  step 2: u(:,n+1)====================================================
k1=aa*tau/(2*hh);
k2=(cc+aa*(hh^2)/6)*tau/(2*hh^3);
k3=(dd-(hh^2)/3)/hh^2;
k4=(ee-dd*(hh^2)/4)/hh^4;
k5=(ff-cc*(hh^2)/12)*tau/(2*hh^5);
k6=ee/(6*hh^4);
k7=bb*tau/(hh*(pp+1));

ld1=k5+k6;
ld2=-k2+k4-4*k5-6*k6;
ld3=-k1+2*k2-k3-4*k4+5*k5+15*k6;
ld4=1+2*k3+6*k4-20*k6;
ld5=k1-2*k2-k3-4*k4-5*k5+15*k6;
ld6=k2+k4+4*k5-6*k6;
ld7=-k5+k6;

A=zeros(J-5,J-5);B=A;
for kk=2:N
    for j=1:J-5
        A(j,j)=ld4;     B(j,j)=ld4;
    end
    for j=1:J-6
        A(j,j+1)=ld5;   B(j,j+1)=ld3;
    end
    for j=1:J-7
        A(j,j+2)=ld6;   B(j,j+2)=ld2;
    end
    for j=1:J-8
        A(j,j+3)=ld7;   B(j,j+3)=ld1;
    end
    for j=2:J-5
        A(j,j-1)=ld3;   B(j,j-1)=ld5;
    end
    for j=3:J-5
        A(j,j-2)=ld2;   B(j,j-2)=ld6;
    end
    for j=4:J-5
        A(j,j-3)=ld1;   B(j,j-3)=ld7;
    end
    for i=1:J-5
        H(i)=0;
        for j=1:J-5
            H(i)=H(i)+B(i,j)*u(j+3,kk-1);
        end
        h(i)=H(i)+k7*(u(i+3-2,kk)^(pp+1)+4*u(i+3-1,kk)^(pp+1)-...
            4*u(i+3+1,kk)^(pp+1)-u(i+3+2,kk)^(pp+1))/6;
    end
    % solve a sevendiaginal system Ax=f
    c=diag(A);
    d=[diag(A,1);0];
    e=[diag(A,2);0;0];
    f=[diag(A,3);0;0;0];
    b=[0;diag(A,-1)];
    a=[0;0;diag(A,-2)];
    g=[0;0;0;diag(A,-3)];
    
    n=length(h);
    
    alpha(1)=c(1);          
    beta(1)=d(1)/alpha(1);
    q(1)=e(1)/alpha(1); 
    k(1)=f(1)/alpha(1);

    gama(2)=b(2);
    alpha(2)=c(2)-gama(2)*beta(1);
    beta(2)=(d(2)-gama(2)*q(1))/alpha(2);
    q(2)=(e(2)-gama(2)*k(1))/alpha(2);
    k(2)=f(2)/alpha(2);

    z(3)=a(3);
    gama(3)=b(3)-z(3)*beta(1);
    alpha(3)=c(3)-z(3)*q(1)-gama(3)*beta(2);
    beta(3)=(d(3)-z(3)*k(1)-gama(3)*q(2))/alpha(3);
    q(3)=(e(3)-gama(3)*k(2))/alpha(3);
    k(3)=f(3)/alpha(3);

    for i=4:n
        m(i)=g(i);
    end

    for i=4:n-3
        z(i)=a(i)-m(i)*beta(i-3);
        gama(i)=b(i)-m(i)*q(i-3)-z(i)*beta(i-2);
        alpha(i)=c(i)-m(i)*k(i-3)-z(i)*q(i-2)-gama(i)*beta(i-1);
        beta(i)=(d(i)-z(i)*k(i-2)-gama(i)*q(i-1))/alpha(i);
        q(i)=(e(i)-gama(i)*k(i-1))/alpha(i);
        k(i)=f(i)/alpha(i);
    end

    z(n-2)=a(n-2)-m(n-2)*beta(n-5);
    gama(n-2)=b(n-2)-m(n-2)*q(n-5)-z(n-2)*beta(n-4);
    alpha(n-2)=c(n-2)-m(n-2)*k(n-5)-z(n-2)*q(n-4)-gama(n-2)*beta(n-3);
    beta(n-2)=(d(n-2)-z(n-2)*k(n-4)-gama(n-2)*q(n-3))/alpha(n-2);
    q(n-2)=(e(n-2)-gama(n-2)*k(n-3))/alpha(n-2);

    z(n-1)=a(n-1)-m(n-1)*beta(n-4);
    gama(n-1)=b(n-1)-m(n-1)*q(n-4)-z(n-1)*beta(n-3);
    alpha(n-1)=c(n-1)-m(n-1)*k(n-4)-z(n-1)*q(n-3)-gama(n-1)*beta(n-2);
    beta(n-1)=(d(n-1)-z(n-1)*k(n-3)-gama(n-1)*q(n-2))/alpha(n-1);

    z(n)=a(n)-m(n)*beta(n-3);
    gama(n)=b(n)-m(n)*q(n-3)-z(n)*beta(n-2);
    alpha(n)=c(n)-m(n)*k(n-3)-z(n)*q(n-2)-gama(n)*beta(n-1);

    % Ly=H
    y(1)=h(1)/alpha(1);
    y(2)=(h(2)-gama(2)*y(1))/alpha(2);
    y(3)=(h(3)-z(3)*y(1)-gama(3)*y(2))/alpha(3);  
    for i=4:n
        y(i)=(h(i)-m(i)*y(i-3)-z(i)*y(i-2)-gama(i)*y(i-1))/alpha(i);
    end

    % Ux=y
    as(n)=y(n);
    as(n-1)=y(n-1)-beta(n-1)*as(n);
    as(n-2)=y(n-2)-beta(n-2)*as(n-1)-q(n-2)*as(n);
    for i=n-3:-1:1
        as(i)=y(i)-beta(i)*as(i+1)-q(i)*as(i+2)-k(i)*as(i+3);
    end
    for j=4:J-2
        u(j,kk+1)=as(j-3);
    end
end
toc;
% %==========error value ============================================
for j=1:J+1
    orer(j)=abs(u(j,N+1)-uu(j,N+1));
end
or2=0;
for j=1:J+1
    or2=or2+(orer(j))^2;
end
or1=norm(orer,inf) % infinite norm
or2=sqrt(hh*or2)   % 2-norm
%========================================================================
%==============mass conservation=====================
Q=0;
for j=2:J
    if N==0
        Q=Q+hh*u(j,1);
    end
    if N>1
        Q=Q+hh/2*(u(j,N)+u(j,N+1));
    end
end
Q
%==============energy conservation=====================
% E1 compute
sum11=0;sum12=0;sum13=0;
for j=2:J
    sum11=sum11+hh*((u(j,2))^2+(u(j,1))^2);
    sum12=sum12+1/hh*((u(j+1,2)-u(j,2))^2+(u(j+1,1)-u(j,1))^2);
    sum13=sum13+1/hh^3*((u(j+1,2)-2*u(j,2)+u(j-1,2))^2+...
        (u(j+1,1)-2*u(j,1)+u(j-1,1))^2);
end
sum14=0;
for j=2:J-1
    sum14=sum14+1/hh^5*((u(j+2,2)-3*u(j+1,2)+3*u(j,2)-u(j-1,2))^2+...
        (u(j+2,1)-3*u(j+1,1)+3*u(j,1)-u(j-1,1))^2);
end
sum14=sum14+1/hh^5*((-3*u(J+1,2)+3*u(J,2)-u(J-1,2))^2+(-3*u(J+1,1)+...
    3*u(J,1)-u(J-1,1))^2);
E1=1/2*sum11+s2/2*sum12+s3/2*sum13-hh^2/12*sum14
% En compute
sum21=0;sum22=0;sum23=0;
for j=2:J
    sum21=sum21+hh*((u(j,N+1))^2+(u(j,N))^2);
    sum22=sum22+1/hh*((u(j+1,N+1)-u(j,N+1))^2+(u(j+1,N)-u(j,N))^2);
    sum23=sum23+1/hh^3*((u(j+1,N+1)-2*u(j,N+1)+u(j-1,N+1))^2+...
        (u(j+1,N)-2*u(j,N)+u(j-1,N))^2);
end
sum24=0;
for j=2:J-1
    sum24=sum24+1/hh^5*((u(j+2,N+1)-3*u(j+1,N+1)+3*u(j,N+1)-...
        u(j-1,N+1))^2+(u(j+2,N)-3*u(j+1,N)+3*u(j,N)-u(j-1,N))^2);
end
sum24=sum24+1/hh^5*((-3*u(J+1,N+1)+3*u(J,N+1)-u(J-1,N+1))^2+...
    (-3*u(J+1,N)+3*u(J,N)-u(J-1,N))^2);
En=1/2*sum21+s2/2*sum22+s3/2*sum23-hh^2/12*sum24;
sum25=0;
for kn=2:N
    for j=2:J
        sum25=sum25+tau*(u(j-1,kn+1)+u(j-1,kn-1)+4*u(j,kn+1)+4*u(j,kn-1)+...
            u(j+1,kn+1)+u(j+1,kn-1))*((u(j+1,kn))^(pp+1)-...
            (u(j-1,kn))^(pp+1));
    end
end
En=En+bb/(12*(pp+1))*sum25

%==================solution figure and error figure========================
% figure;
% figure;
% % hh=0.2; tau=0.1; T=10
% plot(x,u(:,1),'k-',x,u(:,N/2+1),'r--',x,u(:,N+1),'b-.','LineWidth',2);
% for j=1:J+1
%     T5uer(j)=abs(u(j,N/2+1)-uu(j,N/2+1));
%     T10uer(j)=abs(u(j,N+1)-uu(j,N+1));
% end
% figure;
% plot(x,T5uer,'r--',x,T10uer,'b-.','LineWidth',2);

