function progBar = ProgressBar(prog)

FULL_BAR_LENGTH = 20;

numBars = floor(prog*FULL_BAR_LENGTH);
progBar = ['[' 									...
	repmat('=', [1, numBars]) 					...
	repmat(' ', [1, FULL_BAR_LENGTH - numBars]) ...
	']'];

end
