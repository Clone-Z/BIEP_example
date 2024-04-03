%% Polynomial interpolation on a polygon
clear, clc
% Setting the polygon
zv = exp(-1i*pi/4).*[0,1,1+0.5i,0.5+0.5i,0.5+1i,1i,0];

% Setting up the objective functions and test points
f1=@(z) sqrt(z+0.2); f2=@(z) 1./(z.^2+0.04); f3=@(z) 1./(z-1);
zz = chebpts(1e3,[0,1]); 
zz = zz*(zv(2:end)-zv(1:(end-1)))+ones(1e3,6)*diag(zv(1:(end-1))); 
zz = zz(:); F1=f1(zz); F2=f2(zz); F3=f3(zz);

% Discretisation of the boundary of region E
N = 80; a = zv(1:(end-1)); b = zv(2:end);
h = length(a);
ch = chebpts(N,[0,1]); ch = ch.';
t = zeros(h,N);
for k = 1:h
    t(k,:) = (b(k)-a(k))*ch+a(k);
end
t_mid = (t(:,2:end)+t(:,1:(end-1)))./2;
L = abs(t(:,2:end)-t(:,1:(end-1)));

% Solving Symm's equation to obtain the approximate density function
H = ptt_mat(t,t_mid,L,6);
H(:,end)=[]; H(end,:)=[];
F = [zeros(6*(N-1),1);1]; 
U = H\F; U(end)=[];
U = reshape(U,N-1,6);

% Errors at different orders of n
nn = 9:10:299; k=1;
err1 = zeros(size(nn));
err2 = zeros(size(nn));
err3 = zeros(size(nn));
for n=nn
    xi = []; 
    M = round((n+1)*cumsum(diag(L(1:h,:)*U(:,1:h)))); M = [M(1);M(2:end)-M(1:(end-1))];
    for s = 1:h
        temp = den2pts(U(:,s),L(s,:),M(s),a(s),b(s));
        temp = temp(2:end);
        xi = [xi temp]; 
    end
    xi = xi(:);
    wi = BRWeights(xi,[]); wi=wi(:);
    fi1 = f1(xi); fi2 = f2(xi);fi3 = f3(xi);
    err1(k) = norm(bary(zz,fi1,xi,wi)-F1,inf);
    err2(k) = norm(bary(zz,fi2,xi,wi)-F2,inf);
    err3(k) = norm(bary(zz,fi3,xi,wi)-F3,inf); k=k+1;
end

% Plot the error image with respect to order n
figure(1)
semilogy(nn,err1,'ks-', 'linewidth', 1.1, 'markerfacecolor', [56, 109, 225]/255), hold on
semilogy(10:80,exp(-0.418).^(10:80),'--', 'linewidth', 1.1,'Color',[56, 109, 225]/255), hold on
semilogy(nn,err2,'k^-', 'linewidth', 1.1, 'markerfacecolor', [156, 109, 25]/255), hold on
semilogy(10:140,exp(-0.2248).^(10:140),'--', 'linewidth', 1.1,'Color',[156, 109, 25]/255), hold on
semilogy(nn,err3,'ko-', 'linewidth', 1.1, 'markerfacecolor', [156, 9, 225]/255), hold on
semilogy(10:240,exp(-0.1115).^(10:240),'--', 'linewidth', 1.1,'Color',[156, 9, 225]/255), grid on
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')
legend('$f(z)=\sqrt{z+0.2}$','','$f(z)=1/(z^2+0.04)$','','$f(z)=1/(z-1)$','Interpreter','latex')
