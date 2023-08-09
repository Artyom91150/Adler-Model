function F = F2func(Phi, k, Delta, Alpha)

F = 1 ./ (1 + exp(k .* (cos(Delta / 2) - cos(Phi - Alpha - Delta ./ 2))));

end

