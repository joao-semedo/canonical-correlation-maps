function str = EscapeSpecialChars(str)

specialChar	= ['?','%','\'];
rep			= {'''','%%','\\'};

i = 1;
while i <= length(str)
	specialCharIdx = find(str(i) == specialChar);
	if ~isempty(specialCharIdx)
		str = [str(1:i-1) rep{specialCharIdx} str(i+1:end)];
		i = i+length(rep{specialCharIdx});
	else
		i = i+1;
	end
end

end
