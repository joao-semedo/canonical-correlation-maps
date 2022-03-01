function spikes = ZScoreSpikes(spikes, stim)

stimIds	= unique(stim)';
numStim = numel(stimIds);
numPops = numel(spikes);

for stimIdx = 1:numStim
    
    for popIdx = 1:numPops
        spikes{popIdx}(:,:,stim == stimIds(stimIdx)) ...
            = ZScore( spikes{popIdx}(:,:,stim == stimIds(stimIdx)) );
    end
    
end

end



function spikes = ZScore(spikes)

X = Squash(spikes)';
m = mean(X);
s = std(X);
idxs = find( abs(s) < sqrt( eps(class(s)) ) );
if any(idxs)
	s(idxs) = 1;
end

[~, T, numTrials] = size(spikes);
M = repmat( m', [1, T, numTrials] );
S = repmat( s', [1, T, numTrials] );

spikes = (spikes - M) ./ S;

end

