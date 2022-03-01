function firingRates = ComputeFiringRates(spikes, binWidth)

if ~iscell(spikes)
	c{1} = spikes;
	spikes = c;
end

numPops = numel(spikes);
firingRates = cell(1, numPops);
for popIdx = 1:numPops
	firingRates{popIdx} ...
		= (10^3/binWidth) * mean( Squash(spikes{popIdx}), 2 );
end

end
