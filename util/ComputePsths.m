function [psths, resids] = ComputePsths(spikes, stim)

numPops = numel(spikes);

numUnits = zeros(numPops, 1);
for popIdx = 1:numPops
	numUnits(popIdx) = size(spikes{popIdx}, 1);
end

[~, trialLength, numTrials] = size(spikes{1});

stimIds = unique(stim); numStimIds = numel(stimIds);

psths 	= cell(1, numPops);
resids 	= cell(1, numPops);
for popIdx = 1:numPops

	psths{popIdx} 	= zeros(numUnits(popIdx), trialLength, numStimIds);
	resids{popIdx} 	= zeros(numUnits(popIdx), trialLength, numTrials);
	for stimIdx = 1:numStimIds

		psths{popIdx}(:,:,stimIdx) ...
			= mean(spikes{popIdx}(:,:,stim==stimIds(stimIdx)), 3);

		resids{popIdx}(:,:,stim==stimIds(stimIdx)) 			...
			= spikes{popIdx}(:,:,stim==stimIds(stimIdx)) 	...
				- repmat( psths{popIdx}(:,:,stimIdx), 		...
					[1, 1, sum(stim==stimIds(stimIdx))] );

	end

end

end
