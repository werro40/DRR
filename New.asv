%initializing variables
clc
clear
s = 200;%matrix size
t=s/50;
n=1:t:s;
H = zeros(s); %Hamiltonian (normalized with ground state energy)
d = eye(s); %kronecker delta
p = 500; %rho

for i = 1:s
    for j = 1:s
        H(i,j) = d(i,j)*(i^2+p/5+4*p/5*sin(0.8*i*pi)/(0.8*i*pi))-(1-d(i,j))*4*p/5*(sin(0.4*(i-j)*pi)/(0.4*(i-j+d(i,j))*pi)-sin(0.4*(i+j)*pi)/(0.4*(i+j)*pi))*(((-1)^(i+j)+1)/2);
    end
end
%  figure
%  plot(1:s,eig(H),'-r')
% title('Energy Eigenvalues of V_0/E_1=2e6 double well')
% xlabel('n')
% ylabel('$$\frac{E_n}{E_1}$$','interpreter','latex')
x = linspace(0,1,100);
[V,D] = eig(H);
v = V(:,2);
y1 = v(1)*sqrt(2)*sin(pi*x);
for i=2:50
    y1= y1 + v(i)*sqrt(2)*sin(i*pi*x);
end
figure;
plot(x,abs(y1).^2,'r')
title('Probability Density of Second Wavefunction')
xlabel('x/a')
ylabel('$$a^{\frac{1}{2}}\Psi(x)$$','interpreter','latex')