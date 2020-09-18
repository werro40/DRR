x%initializing variables
clc
clear
s = 200;%matrix size
t=s/50;
n=1:t:s;
H = zeros(s); %Hamiltonian (normalized with ground state energy)
d = eye(s); %kronecker delta
p = 50; %rho

for i = 1:s
    for j = 1:s
        H(i,j) = d(i,j)*(j^2+(pi*p)^2/48*(1-6/(pi*j)^2))+(1-d(i,j))*((p^2/4)*((-1)^(i+j)+1)*(1/(i-j+d(i,j))^2-1/(i+j)^2));
    end
end
% figure
% plot(1:s,eig(H),'-r',n,50*n-25,'-k',n,1*n.^2,'--g',n,n.^2+(50*pi)^2/48,':bs')
% title('Energy Eigenvalues of Harmonic Oscillator')
% xlabel('n')
% ylabel('$$\frac{E_n}{E_1}$$','interpreter','latex')
x = linspace(0,1,100);
[V,D] = eig(H);
v = V(:,2);
y1 = v(1)*sqrt(2)*sin(pi*x);
y2 = v(1)*sqrt(2)*sin(pi*x);
y3 = v(1)*sqrt(2)*sin(pi*x);
for i=2:10
    y1= y1 + v(i)*sqrt(2)*sin(i*pi*x);
end
for i=2:20
    y2= y2 + v(i)*sqrt(2)*sin(i*pi*x);
end
for i=2:50
    y3= y3 + v(i)*sqrt(2)*sin(i*pi*x);
end
y= 5*sqrt(2)*pi*(25*pi)^0.25*exp(-25/2*pi^2*(x-0.5).^2).*(1/2-x);
figure;
plot(x,y1,'r',x,y2,'*-g',x,y3,'--b',x,y,'-dk')
title('First Excited State Wavefunction of Harmonic Oscillator')
xlabel('x/a')
ylabel('$$a^{\frac{3}{2}}\Psi(x)$$','interpreter','latex')
