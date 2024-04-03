function H = ptt_mat(Zi,Zi_mid,Li,N)
%m个区域相互作用
%N个E
[m,~] = size(Zi);
T = cell(m);
for i=1:m
    for j=1:m
        T{i,j}=(-1)^((j<=N)+1)*p_matrix(Zi_mid(i,:),Zi(j,:),Zi_mid(j,:),Li(j,:));
    end
end
H = cell2mat(T);
LE = Li(1:N,:); LF = Li((N+1):end,:);
LE = LE.'; LF = LF.';
LE = LE(:); LF=LF(:);
H1 = [-ones(size(LE)),zeros(size(LE));zeros(size(LF)),ones(size(LF))];
H2 = [LE,zeros(size(LE));zeros(size(LF)),LF]; H2=H2.';
H = [H,H1;H2,zeros(2,2)];
end

function h = p_matrix(x_mid,z,z_mid,li)
%x_mid 测试点
%z 线段端点
%z_mid 线段中点
%li线段长
x_mid = x_mid(:);
z = z(:); z = z.';
z_mid = z_mid(:); z_mid = z_mid.';
li = li(:); li = li.';
if prod(x_mid.'==z_mid)
h = -(log(abs(bsxfun(@minus,x_mid,z(1:(end-1))))) ...
    +4*log(abs(bsxfun(@minus,x_mid,z_mid))+eye(length(li))) ...
    +log(abs(bsxfun(@minus,x_mid,z(2:end)))))*diag(li)/6;

h(logical(eye(length(li)))) = li.*(1-log(li/2));
else
    h = -(log(abs(bsxfun(@minus,x_mid,z(1:(end-1))))) ...
    +4*log(abs(bsxfun(@minus,x_mid,z_mid))) ...
    +log(abs(bsxfun(@minus,x_mid,z(2:end)))))*diag(li)/6;
end
end