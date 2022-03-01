function numWorkers = StartParPool(numParRuns)

currParPool = gcp('NoCreate');

numWorkers = min(feature('NumCores'), numParRuns);

if numWorkers > 0
    if isempty(currParPool)
        parpool( numWorkers );
    elseif currParPool.NumWorkers < numWorkers
        delete(currParPool)
        parpool( numWorkers );
    end
end

end
