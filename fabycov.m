function fabycov(X)


X = [];

error(nargchk(1,1,nargin))

[r c] = size(X);

if r ~= c,
    
    error('Input matrix is not square.');
    
elseif ~all(all(X == X')),
    
    error('Input matrix is not symmetric.');
    
elseif any(diag(X) <= 0),
    
    error('The covariance/correlation matrix must be positive definite.');
    
else
    
end

dX = diag(X);
D = 1./sqrt(dX);
R = X.*(D*D');
p = r;
[A L A] = svd(R,0);
l = diag(L);
P = (l / sum(l))*100;
CP = cumsum(P);
ft = [1:p]';

disp(' ')
disp('Table of General Structure of Extracted Components.')
fprintf('-------------------------------------------------------------\n');
disp('                            Percent of       Cummulative       ');
disp(' Factors     Eigenvalue      Variance     Percent of Variance  '); 
fprintf('-------------------------------------------------------------\n');
fprintf('    %d      %10.4f      %10.4f        %10.4f\n',[ft,l,P,CP].');
fprintf('-------------------------------------------------------------\n');

B = A*sqrt(L);
f = L >= 1.0;
m = sum(sum(f));

disp('  ');
fprintf('By benchmarking the latent root citerion, the number of retainable factors suggested are %.i\n', m);
disp('  ');
ask = input('Do you need to work with this number of factors? (y/n): ','s');

if~strcmp(ask,'y'),
    
   m = input('Enter the number of factors you need: ');
   
   while (m > p)
       
       disp(' ');
       fprintf('Error: The number of factors requested is too large for the number of the observed variables. It must be equal or lesser than %.i\n', p);
       disp(' ');
       m = input('Enter the number of factors you need: ');
       
   end
   
end

F = B(:,1:m);
pt = sum(F.^2)/p;
pp = [pt sum(pt)];

C = F.^2;
C = sum(C,2);
Factors = [F C];

disp(' ')
disp('Table of Unrotated Principal Components of the Factor Analysis.')
disp('-----------------------------------------------------------------------------');
Factors
fprintf('-----------------------------------------------------------------------------\n');
fprintf('On factors, Factor 1 = column 1 and so forth to %.i\n', m);
disp('Values on the last column are the Communality');
fprintf('On variates, Variate 1 = first row and so forth to %.i\n', p);
disp(' ')
disp('Proportion of Total (Standardized) Sample Variance.');
disp('-----------------------------------------------------------------------------');
pp
fprintf('-----------------------------------------------------------------------------\n');

sv = 1 - C;
Re = F*F' + diag(sv);
Rm = R - Re;

disp('  ');
rt = input('Do you want to apply a varimax factor rotation? (y/n): ','s');

if rt == 'y'
    
    loadings = F;
    b = loadings;
    [n,nf] = size(loadings);
    hjsq = diag(loadings*loadings');
    hj = sqrt(hjsq);

    for iter = 1:10
        
        for i = 1:nf-1,
            
            jl = i + 1;
            
            for j = jl:nf,
                
                xj = loadings(:,i)./hj;
                yj = loadings(:,j)./hj;
                uj = (xj.*xj) - (yj.*yj);
                vj = 2*xj.*yj;
                A = sum(uj);
                B = sum(vj);
                C = (uj'*uj) - (vj'*vj);
                D = 2*uj'*vj;
                num = D - 2*A*B/n;
                den = C - (A^2 - B^2)/n;
                tan4p = num/den;
                phi = atan2(num,den)/4;
                angle = phi*180/pi;
                
                if abs(phi) > 0.00001;
                    
                    Xj = (cos(phi)*xj) + (sin(phi)*yj);
                    Yj = (-sin(phi)*xj) + (cos(phi)*yj);
                    bj1 = Xj.*hj;
                    bj2 = Yj.*hj;
                    b(:,i) = bj1;
                    b(:,j) = bj2;
                    loadings(:,i) = b(:,i);
                    loadings(:,j) = b(:,j);
                    
                end
                
            end
            
        end
        
        loadings = b;
        
    end
    
    F = loadings;
    pt = sum(F.^2)/p;
    pp = [pt sum(pt)];
    
    C = F.^2;
    C = sum(C,2);
    Factors = [F C];
    
    disp(' ')
    disp('Table of Varimax Rotated Principal Components of the Factor Analysis.')
    disp('-----------------------------------------------------------------------------');
    Factors
    fprintf('-----------------------------------------------------------------------------\n');
    fprintf('On factors, Factor 1 = column 1 and so forth to %.i\n', m);
    disp('Values on the last column are the Communality');
    fprintf('On variates, Variate 1 = first row and so forth to %.i\n', p);
    disp(' ')
    disp('Table of Cumulative Proportion of Total (Standardized) Sample Variance.');
    disp('-----------------------------------------------------------------------------');
    pp
    fprintf('-----------------------------------------------------------------------------\n');
    
    sv = 1 - C;
    Re = (F*F') + (diag(sv));
    Rm = R - Re;
    
else
    
end

disp(' ');
rm = input('Do you need to output the residual matrix? (y/n): ','s');
disp(' ');

if rm == 'y',
    
    disp('Residual matrix:');
    Rm
    
else
    
end

return
