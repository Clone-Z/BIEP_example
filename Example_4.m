%% Rational interpolation of branch singularities
clear,clc
% Setting up the objective function and test points
f = @(x) exp(1./sqrt(1e6*(x).^2+1));
zz = linspace(-1,1,1e5);zz = zz(:); ff=f(zz);

% Pre-solve the discrete points t_E of E
a = [-1;0.001i;-0.001i]; b = [1;10.001i;-10.001i]; N = 300;
X = chebpts(N,[0,1]); X = X.';
t = a+(b-a)*X;
t_mid = (t(:,1:(end-1))+t(:,2:end))/2;
L = abs(-t(:,1:(end-1))+t(:,2:end));
H = ptt_mat(t,t_mid,L,1); 
F = [zeros(3*(N-1),1);1;1]; 
U = H\F; U(end)=[]; U(end)=[];
U = reshape(U,N-1,3); 
t(1,:) = den2pts(U(:,1),L(1,:),N-1,a(1),b(1));
t(2,:) = den2pts(U(:,2),L(2,:),N-1,a(2),b(2));
t(3,:) = den2pts(U(:,2),L(2,:),N-1,a(3),b(3));

% Solving Symm's equation to obtain the approximate density function
t_mid = (t(:,1:(end-1))+t(:,2:end))/2;
L = abs(-t(:,1:(end-1))+t(:,2:end));
H = ptt_mat(t,t_mid,L,1); 
U = H\F; c1 = U(end); c2 = U(end-1);
U(end)=[]; U(end)=[]; U = reshape(U,N-1,3); 

% Errors at different orders of n
nn = 4:4:200; err = zeros(size(nn)); k=1;
for n = nn
    M(1) = floor(n/2); M(2) = n-M(1);
    xi = den2pts(U(:,1),L(1,:),n,a(1),b(1)); xi = xi(:);
    pj1 = den2pts(U(:,2),L(2,:),M(1)-1,a(2),b(2));
    pj2 = den2pts(U(:,3),L(3,:),M(2)-1,a(3),b(3));
    pj = [pj1,pj2];pj = pj(:);
    wi = BRWeights(xi,pj); 
    fi = f(xi);
    err(k) = norm((bary(zz,fi,xi,wi)-ff),inf); k=k+1;
end

semilogy(nn,err,'ko-', 'linewidth', 1.1, 'markerfacecolor', [156, 9, 225]/255), hold on
semilogy(nn,1e2*exp(-(c1+c2)*nn),'--', 'linewidth', 1.1), grid on
axis([0 200 1e-16 1e0])
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')
