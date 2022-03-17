function r = CanonCorrFitAndPredict(Xtrain, Ytrain, Xtest, Ytest)
warning('off','stats:canoncorr:NotFullRank');
lastwarn('')

[A, B, ~] = canoncorr( Xtrain, Ytrain );
[~, id] = lastwarn;

numPairs = min( size(Xtrain, 2), size(Ytrain, 2) );
if ~isempty(id) && strcmp(id, 'stats:canoncorr:NotFullRank')
	r = nan(1, numPairs);
else
	r = zeros(1, numPairs);
	for pairIdx = 1:numPairs
		r(pairIdx) = corr( Xtest*A(:,pairIdx), Ytest*B(:,pairIdx) );
	end
end

warning('on','stats:canoncorr:NotFullRank');
end
