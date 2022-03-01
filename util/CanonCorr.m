function [A, B, r, warn] = CanonCorr(X,Y)
warning('off','stats:canoncorr:NotFullRank');
lastwarn('')

[A, B, r] = canoncorr( X, Y );
[~, id] = lastwarn;

warn = false;
if strcmp(id, 'stats:canoncorr:NotFullRank')
    warn = true;
end

warning('on','stats:canoncorr:NotFullRank');
end
