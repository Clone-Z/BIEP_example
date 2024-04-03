function xi = den2pts(u,li,M,a,b)
% u:the numble of density
% li:the length of interval
% M:the quantity of points
% [a,b]:the domain of points
%M = M-1;
li = li(:); li = li.';
if M==0
    xi = [a,b];
else
u = u(~(sign(u)+1==0)); u = u/sum(u);
li = li(~(sign(u)+1==0));
u = u(:); u = u';
% mu = sum(tril(repmat(u,length(u),1)),2);
mu = [0,cumsum(u.*li)];mu = mu(:);
t = mu(end)/M*(0:M);
A = (mu(1:(end-1))-t)./(repmat(mu(1:(end-1))-mu(2:end),1,length(t)));
B = (sign(A-1)+1)/2; C = (sign(A)+1)/2;
D = diag(li)*(A.*(1-B).*C+B);
xi = sum(D)/sum(D(:,end));
xi = (b-a)*xi+a;
end

% uk = prod( 1./(xi'-xi+eye(M+1)),2);
% uk = uk/max(uk);

end