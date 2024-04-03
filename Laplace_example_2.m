%% Rational interpolation to solve the Laplace equation on a smooth boundary region
clear, clf
% Boundary condition and test points
f = @(x) (sin(3*imag(log(x))));
G = @(the) (1.5+0.2*cos(5*the)).*exp(1i*the);
zz = G(linspace(0,2,1e3)*pi); zz = zz(:);
ff = f(zz);

% The discretisation of the boundary of E
N = 300; the = linspace(0,2,N)*pi; te = G(the); 
[the_mid,Lthe]  = Segm(the); te_mid = G(the_mid); 
Le = abs(te(1:(end-1))-te_mid)+abs(te(2:end)-te_mid);

% Setting up the F-region
load five_poles.mat                      % First 5 poles calculated by AAA
b = (abs(a)+2)./abs(a).*a;
z = chebpts(N,[0,1]); 
tf = zeros(length(a),N); 
for s = 1:length(a)
    tf(s,:) = (b(s)-a(s))*z+a(s);
end
[tf_mid,Lf]  = Segm(tf);
t = [te;tf]; L = [Le;Lf];
t_mid = [te_mid;tf_mid];

% Solving Symm's equation to obtain the approximate density function
H = ptt_mat(t,t_mid,L,1); 
F = [zeros(6*(N-1),1);1;1]; 
U = H\F; c2 = U(end); c1=U(end-1);
U(end)=[]; U(end)=[];
U = reshape(U,N-1,6); 
U(:,1) = U(:,1).*(L(1,:).'./Lthe.');
U = U/(Lthe(1,:)*U(:,1));

% Errors at different orders of n
nn = 10:10:300; err2 = zeros(size(nn)); k = 1;
for n = nn
    The = den2pts(U(:,1),Lthe,n+1,0,2)*pi; The = The(2:end);
    xi = G(The); xi = xi(:);

    yThe = [1/3*The(1:(end-1))+2/3*The(2:end) 2/3*The(1:(end-1))+1/3*The(2:end)];
    yThe = yThe(:); yi = G(yThe);

    zi = []; m = diag(L*U); m(1)=[];
    m = round(n*cumsum(m)); m = [m(1);m(2:end)-m(1:(end-1))];
    for s = 1:length(a)
        temp = den2pts(U(:,s+1),L(s+1,:),m(s)-1,a(s),b(s));
        zi = [zi, temp];       % Poles of rational interpolation
    end
    zi = zi(:);
    wk = BRWeights(xi,zi);

    E = eye(length(xi));
    A = real(bary(yi,E,xi,wk));
    B = imag(bary(yi,E,xi,wk)); 
    F = (A*f(xi))-f(yi); 
    [S,C,V] = svd(B,'econ');
    aj = V*((S'*F)./diag(C));    % Approximation of the imaginary part
    
    fz = f(xi)+1i*aj;
    R = real(bary(zz,fz,xi,wk)); 
    err2(k) = norm((ff-R),inf); k = k+1;
end

figure(1)
plot(real(xi),imag(xi),'.'), hold on
plot(real(zi),imag(zi),'.'),axis equal
figure(2)
semilogy(nn,err2,'ko-', 'linewidth', 1.1, 'markerfacecolor', [156, 9, 225]/255),hold on
semilogy(nn,1e2*exp(-nn*(c1+c2)),'--r', 'linewidth', 1.1), grid on
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')
axis([0 300 1e-16 1e0])

function [Zi_mid,Li] = Segm(Zi)
Zi_mid = (Zi(:,2:end)+Zi(:,1:(end-1)))./2;
Li = abs(Zi(:,2:end)-Zi(:,1:(end-1)));
end