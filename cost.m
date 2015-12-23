function F = cost(y, yStar)
	F = sum((y - yStar) .^ 2);
end
