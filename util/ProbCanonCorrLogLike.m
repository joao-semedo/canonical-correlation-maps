function logLike = ProbCanonCorrLogLike(X, Y, estParams)

numelRank = numel(estParams);
logLike = zeros( 1, numelRank );
for i = 1:numelRank;
    W1		= estParams(i).W1;
    W2		= estParams(i).W2;
    Psi1	= estParams(i).Psi1;
    Psi2	= estParams(i).Psi2;
    mu1		= estParams(i).mu1;
    mu2		= estParams(i).mu2;
    
    S = [W1*W1'+Psi1,W1*W2';W2*W1',W2*W2'+Psi2];
    m = [mu1; mu2]';
    
    logLike(i) = MvnLogLike([X Y], m, S);
end

end

