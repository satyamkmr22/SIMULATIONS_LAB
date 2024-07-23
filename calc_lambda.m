function lambda = calc_lambd(Bi)
    initial_guess = linspace(1, 200, 200);
    result = zeros(200, 1);
    
    for j = 1:numel(initial_guess)
        fun = @(x) (x * tan(x)) - Bi;
        x = fsolve(fun, initial_guess(j));
        result(j) = x;
    end
    lambda=result;
end
