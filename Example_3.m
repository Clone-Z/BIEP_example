%% Rational interpolation of isolated singularities
clear,clc
% Setting up the objective function and test points
f = @(x) exp(1./(1e4*(x).^2+1));
zz = linspace(-1,1,1e5); zz = zz(:); ff=f(zz);

% Pre-solve the discrete points t_E of E
a = -1; b = 1; N = 100;
te = a+(b-a)*chebpts(N,[0,1]); te = te.';
rho = linspace(-pi/2,3*pi/2,N);
rho_mid = (rho(1:(end-1))+rho(2:end))/2;
tf = 1e-3*exp(1i*rho)+0.01i; 
tf_mid = 1e-3*exp(1i*rho_mid)+0.01i;
t = [te;tf;(tf').'];
t_mid = [(te(:,1:(end-1))+te(:,2:end))/2;tf_mid;(tf_mid').'];
L = [abs(-te(:,1:(end-1))+te(:,2:end));0.2*pi/(N-1)*ones(1,N-1);0.2*pi/(N-1)*ones(1,N-1)];
H = ptt_mat(t,t_mid,L,1); 
F = [zeros(3*(N-1),1);1;1]; 
U = H\F; U(end)=[]; U(end)=[];
U = reshape(U,N-1,3); 
te = den2pts(U(:,1),L(1,:),N-1,-1,1);

% Solving Symm's equation to obtain the approximate density function
t = [te;tf;(tf').'];
t_mid = [(te(:,1:(end-1))+te(:,2:end))/2;tf_mid;(tf_mid').'];
L = [abs(-te(:,1:(end-1))+te(:,2:end));0.2*pi/(N-1)*ones(1,N-1);0.2*pi/(N-1)*ones(1,N-1)];
H = ptt_mat(t,t_mid,L,1); 
F = [zeros(3*(N-1),1);1;1]; 
U = H\F; c1 = U(end); c2 = U(end-1);
U(end)=[]; U(end)=[];
U = reshape(U,N-1,3); 

% Errors at different orders of n
nn = 4:4:120; err = zeros(size(nn)); k=1;
for n = nn
    M = round(n*cumsum(diag(L(2:3,:)*U(:,2:3)))); M = [M(1);M(2:end)-M(1:(end-1))];
    xi = den2pts(U(:,1),L(1,:),n,-1,1); xi = xi(:);
    pj1 = den2pts(U(:,2),L(2,:),M(1),-pi/2,3*pi/2); pj1(end)=[];
    pj2 = den2pts(U(:,3),L(3,:),M(2),pi/2,-3*pi/2); pj2(end)=[];
    pj = [1e-3*exp(1i*pj1)+0.01i,1e-3*exp(1i*pj2)-0.01i];pj = pj(:);
    wi = BRWeights(xi,pj); 
    fi = f(xi);
    err(k) = norm((bary(zz,fi,xi,wi)-ff),inf); k=k+1;
end

semilogy(nn,err,'ko-', 'linewidth', 1.1, 'markerfacecolor', [156, 9, 225]/255), hold on
semilogy(nn,1e3*exp(-(c1+c2)*nn),'--', 'linewidth', 1.1), grid on
axis([0 120 1e-16 1e0])
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')

