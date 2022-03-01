function spikes = BinTime(spikeRasters, binWidth)

spikes = spikeRasters;

if binWidth == 1
	return
end

[numUnits, trialLength, numTrials] = size(spikeRasters);

binnedTrialLength = floor(trialLength/binWidth);

if mod(trialLength, binWidth) ~= 0
	spikes(:,binnedTrialLength*binWidth + 1:end,:) = [];
end

spikes = permute(spikes, [2 1 3]);
spikes = reshape(spikes, ...
	[binWidth, numUnits*binnedTrialLength, numTrials]);
spikes = sum(spikes, 1);
spikes = reshape(spikes, [binnedTrialLength, numUnits, numTrials]);
spikes = permute(spikes, [2 1 3]);

end
