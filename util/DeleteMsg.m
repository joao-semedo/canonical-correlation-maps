function [] = DeleteMsg(msg)

msgLength = length(msg);
if msgLength > 0
	fprintf( repmat('\b', [1, msgLength]) );
end

end
