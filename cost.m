function F = cost(y_ampl, y_per, y)
    F_ampl = sum((y_ampl(floor(length(y_ampl) / 4) : end) - y(1)) .^ 2);
    F_per = sum((y_per(floor(length(y_per) / 4) : end) - y(2)) .^ 2);
    F_povp1 = sum((0.5 - min(abs(mean(y(1)) - y_ampl), 0.5)) .^ 3);
    F_povp2 = sum((0.5 - min(abs(mean(y(2)) - y_per), 0.5)) .^ 3);

    F = F_ampl + F_per + F_povp1 + F_povp2;
end
