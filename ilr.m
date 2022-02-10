function y = ilr(x,Q,type)

[~,n] = size(x);

if nargin < 2
    
    V = eye(n);
    V = V(1:n,1:n-1) - [zeros(1,n-1);V(2:n,2:n)];
    U = local_gramschmidt(V);
    
else
    
    if nargin == 2
        
        type = 'partition';
        
    end
    
    if strcmp(type,'partition')
        
        pos = Q==1;
        neg = Q==-1;
        r = repmat(sum(pos,2),1,n);
        s = repmat(sum(neg,2),1,n);
        U = (pos.*(sqrt(s./(r.*(r + s)))) - neg.*(sqrt(r./(s.*(r + s)))))';
        
    elseif strcmp(type,'basis')
        
        U = Q;
        
    end
    
    if round(norm(U).*1000000) ~= 1000000
        
        warning('ilr:IllDefinedBasis','ilr(): Orthonomal basis is not well-defined using default basis')
        V = eye(n);
        V = V(1:n,1:n-1) - [zeros(1,n-1);V(2:n,2:n)];
        U = local_gramschmidt(V);
        
    end
    
end

y = local_clr(x)*U;

end

function y = local_clr(x)

[~, n] = size(x);
y = log(x./repmat(geomean(x,2),1,n));

end

function [Q,R] = local_gramschmidt(A)

[~,n] = size(A);
R = zeros(n-1,n);
Q = zeros(n+1,n);

for j = 1:n
    
   v = A(:,j);
   
   for i = 1:j-1
       
        R(i,j) = Q(:,i)'*A(:,j);
        v = v - R(i,j)*Q(:,i);
        
   end
   
   R(j,j) = norm(v);
   Q(:,j) = v/R(j,j);
   
end

end
