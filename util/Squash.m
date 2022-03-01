function B = Squash(A)

dimA = size(A);

if numel(dimA) == 3
    B = reshape(A, [dimA(1), dimA(2)*dimA(3)]);
else
    B = A;
end

end
