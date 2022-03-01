function jSpikes = JitterSpikes(spikes, jitterWindow)

[numUnits, T, ~] = size(spikes);

jSpikes = zeros( size(spikes) );

numWindows = ceil(T/jitterWindow);
for windowIdx = 1:numWindows
    for unitIdx = 1:numUnits
        
        unitSpikes = permute( spikes(unitIdx, ...
            (windowIdx-1)*jitterWindow + 1:...
            min( windowIdx*jitterWindow, T ),: ), ...
            [3 2 1]);
        
        [trialIdxs, spikeTimes] = find(unitSpikes > 0);
        
        numSpikes = numel(trialIdxs);
        randIdxs = randperm(numSpikes);
        jTrialIdxs = trialIdxs( randIdxs );
        jSpikeTimes = spikeTimes( randIdxs );
        
        for spikeIdx = 1:numSpikes
            i = unitIdx;
            j = jSpikeTimes(spikeIdx) + (windowIdx-1)*jitterWindow;
            k = trialIdxs(spikeIdx);
            
            jSpikes( i,j,k ) = jSpikes( i,j,k ) ...
                + unitSpikes( jTrialIdxs(spikeIdx), ...
                        jSpikeTimes(spikeIdx) );
        end
        
    end
end

end

