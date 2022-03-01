function estParams = ComputeProbCanonCorrParams(X1, X2, A, B, r, q)

X = [ X1, X2 ];
mu = mean(X, 1)';
S = cov(X);

m1 = size(X1, 2);

numDims = numel(q);
estParams(numDims,1).W1 = 0;
for i = 1:numDims
    
    Rd = diag(r); Rd = Rd(1:q(i),1:q(i));
    
    M1 = Rd.^(1/2); M2 = M1;
    
    W1 = S(1:m1,1:m1)*A(:,1:q(i))*M1;
    W2 = S(m1+1:end,m1+1:end)*B(:,1:q(i))*M2;
    Psi1 = S(1:m1,1:m1)-W1*W1';
    Psi2 = S(m1+1:end,m1+1:end)-W2*W2';
    mu1 = mu(1:m1);
    mu2 = mu(m1+1:end);
    
    estParams(i).W1 = W1;
    estParams(i).W2 = W2;
    estParams(i).Psi1 = Psi1;
    estParams(i).Psi2 = Psi2;
    estParams(i).mu1 = mu1;
    estParams(i).mu2 = mu2;
    
    estParams(i).A = A;
    estParams(i).B = B;
    
end

end
