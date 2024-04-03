%% Polynomial interpolation on disjoint regions
clear,clc
% Discretisation of the boundary of region E
N = 100; 
zv = [-1-0.1i,-1+0.1i,1+0.1i,1-0.1i,-1-0.1i];
zv = zv(:); A = zv(1:(end-1)); B = zv(2:end);
X = chebpts(N,[0,1]); X = X.';
t1 = A+(B-A)*X;
rho = linspace(-pi/2,3*pi/2,N);
rho_mid = (rho(1:(end-1))+rho(2:end))/2;
t2 = 0.15*exp(1i*rho)+0.8i; 
t2_mid = 0.15*exp(1i*rho_mid)+0.8i;
t = [t1;t2;(t2').'];
t_mid = [(t1(:,1:(end-1))+t1(:,2:end))/2;t2_mid;(t2_mid').'];
L = [abs(-t1(:,1:(end-1))+t1(:,2:end));0.3*pi/(N-1)*ones(1,N-1);0.3*pi/(N-1)*ones(1,N-1)];

% Solving Symm's equation to obtain the approximate density function
H = ptt_mat(t,t_mid,L,6); H(:,end)=[]; H(end,:)=[];
F = [zeros(6*(N-1),1);1]; 
U = H\F; Ve = U(end); U(end)=[];
U = reshape(U,N-1,6); 

% Setting up the objective functions and test points
a = 0.5; f = @(x) 1./((x-a).^2+0.2);
zz = chebpts(1e3,[0,1]); zv = zv.';
zz = zz*(zv(2:end)-zv(1:(end-1)))+ones(1e3,4)*diag(zv(1:(end-1))); 
zz = zz(:); zz1 = 0.15*exp(1i*linspace(-pi/2,3*pi/2,1e3))+0.8i; 
zz1 = zz1(:); zz = [zz;zz1;(zz1.')'];
ff=f(zz); 

% Errors at different orders of n
nn = 20:20:1000; err = zeros(size(nn)); k=1;
for n = nn+1
    M = round(n*cumsum(diag(L*U))); M = [M(1);M(2:end)-M(1:(end-1))];
    xi1 = den2pts(U(:,1),L(1,:),M(1),A(1),B(1)); xi1(end)=[];
    xi2 = den2pts(U(:,2),L(2,:),M(2),A(2),B(2)); xi2(end)=[];
    xi3 = den2pts(U(:,3),L(3,:),M(3),A(3),B(3)); xi3(end)=[];
    xi4 = den2pts(U(:,4),L(4,:),M(4),A(4),B(4)); xi4(end)=[];
    xi5 = 0.15*exp(1i*den2pts(U(:,5),L(5,:),M(5),-pi/2,3*pi/2))+0.8i; xi5(end)=[];
    xi6 = 0.15*exp(1i*den2pts(U(:,6),L(6,:),M(6),pi/2,-3*pi/2))-0.8i; xi6(end)=[];
    xi = [xi1,xi2,xi3,xi4,xi5,xi6];xi = xi(:);
    wi = BRWeights(xi,[]); wi=wi(:);
    fi = f(xi);
    err(k) = norm((bary(zz,fi,xi,wi)-ff),inf); k=k+1;
end
%%
figure(1)
N = length(xi);
g = @(s) prod(s-xi);
u = @(s) -log(abs(g(s)))/(N);
    
xgrid = -1.5:.05:1.5; ygrid = -1.5:.05:1.5;
[xx,yy] = meshgrid(xgrid,ygrid); ss = xx+1i*yy; ss = ss(:);
uss = u(ss.'); uss = reshape(uss,length(ygrid),length(xgrid));
contour(xx,yy,uss,40,'linewidth', 1);
colorbar,hold on
plot(real(xi),imag(xi),'b.','markersize',7), hold on
plot(a,sqrt(0.2),'+r','LineWidth',2,'markersize',7)

figure(2)
c = -Ve+u(a+1i*sqrt(0.2));
semilogy(nn,err,'ko-', 'linewidth', 1.1, 'markerfacecolor', [156, 9, 225]/255), hold on
semilogy(nn,1e4*exp(c*nn),'--', 'linewidth', 1.1), grid on
axis([0 1000 1e-16 1e2])
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')
