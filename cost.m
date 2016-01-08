function F = cost(y, yStar)
    F = sum((y(floor(length(y) / 4) : end) - yStar) .^ 2);
end
