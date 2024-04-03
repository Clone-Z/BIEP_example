function wi = BRWeights(xi,pj)
    xi = xi(:); pj = pj(:);
    D = bsxfun(@minus,xi,xi.')+eye(length(xi)); 
    lgD = round(log(abs(D))); D = D./(exp(lgD));
    N = bsxfun(@minus,xi,pj.');
    lgN = round(log(abs(N))); N = N./(exp(lgN));
    wi = (prod(N,2)./prod(D,2)).*(exp(sum(lgN,2)-sum(lgD,2)));
    wi = wi/norm(wi,inf);
end
