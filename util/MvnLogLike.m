function logLike = MvnLogLike(X, m, S)

[n, p] = size(X);

M = m( ones(n, 1), : );
X = (X - M);

logLike = -(1/2)*( n*p*log(2*pi) + n*logdet(S) ...
	+ sum(sum( X .* ( S\X' )' )) );

end
