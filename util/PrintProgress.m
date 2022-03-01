function msg = PrintProgress(iter, numIter, iterTime)

COMMAND_WINDOW_WIDTH = 80; % columns

prog = iter/numIter;

msg = sprintf([ProgressBar(prog),' ','%3i%%'], floor(prog*100));

hasIterTime = (nargin > 2);
if hasIterTime
	if prog < 1
		eta = mean( iterTime(1:iter) )*(numIter-iter);
		msg = [msg,'     ','ETA: ' ...
			datestr(eta/(24*60*60), 'DD-HH:MM:SS')];
	elseif prog == 1
		elapsedTime = sum(iterTime);
		msg = [msg,'     ','Elapsed Time: ' ...
			datestr(elapsedTime/(24*60*60), 'DD-HH:MM:SS')];
	end
end

msg = [msg repmat( ' ', [1, COMMAND_WINDOW_WIDTH - length(msg)] )];

fprintf(EscapeSpecialChars(msg));

end
