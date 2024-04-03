%% Approximating the harmonic function with the real part of polynomial interpolation
clear;clf
% Setting the polygon
zv = exp(-1i*pi/4).*[0,1,1+0.5i,0.5+0.5i,0.5+1i,1i,0];

% Setting up the objective functions and test points 
f=@(z) log(abs(z-1)); 
zz = chebpts(1e3,[0,1]); 
zz = zz*(zv(2:end)-zv(1:(end-1)))+ones(1e3,6)*diag(zv(1:(end-1))); 
zz = zz(:); ff=f(zz); 

% Discretisation of the boundary of region E
N = 80; a = zv(1:(end-1)); b = zv(2:end);
h = length(a);
ch = chebpts(N,[0,1]); ch = ch.';
t = zeros(h,N);
for s = 1:h
    t(s,:) = (b(s)-a(s))*ch+a(s);
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
nn = 10:10:300; err = zeros(size(nn)); k=1;
for n=nn
    xi = []; 
    M = round((n)*cumsum(diag(L(1:h,:)*U(:,1:h)))); M = [M(1);M(2:end)-M(1:(end-1))];
    for s = 1:h
        temp = den2pts(U(:,s),L(s,:),M(s),a(s),b(s));
        temp = temp(2:end);
        xi = [xi temp]; 
    end
    xi = xi(:);
    xj = [xi(2:end);xi(1)];
    wi = BRWeights(xi,[]); 

    yi = [xi/3+2*xj/3; 2*xi/3+1*xj/3]; 
    E = eye(length(xi));
    A = real(bary(yi,E,xi,wi));
    B = imag(bary(yi,E,xi,wi)); 
    F = (A*f(xi))-f(yi); 
    [S,C,V] = svd(B,'econ');
    aj = V*((S'*F)./diag(C));        % Approximation of the imaginary part

    fz = f(xi)+1i*aj;
    R = real(bary(zz,fz,xi,wi)); 
    err(k) = norm((ff-R),inf); k = k+1;
end
semilogy(nn,err,'ko-', 'linewidth', 1.1, 'markerfacecolor', [156, 9, 225]/255), grid on
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')