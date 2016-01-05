function find_model
    % number of iterations
    max_iterations = 1000;
    % population size
    pop_size = 20;
    
    % create initial population
    % creates cell array of empty matrices
    pop_array = cell(1, pop_size); %cell array
    % protein concentrations
    p_con = cell(1, pop_size);
    for i = 1 : pop_size
        % zacetno stevilo proteinov
        M = ones(3,10);
        M(:, 1) = 0; % spremeni tip izrazanja na gensko izrazanje
        M(:, 5) = 0;
        M(:, 6) = 0;
        M(:, 7) = 0; % spremeni tip degradacije na linearno alfa = 1
        pop_array{i} = M;
        
        p_con{i} = zeros(1, 3);
        % TODO: Zacetne vrednosti.
    end
    
    % create learning examples
    % sinus s frekvenco 10 Hz s 1000 vzorci.
    x = sin(10 * 0 : 0.01 : 9.99);
    y_ampl = zeros(size(x));
    for i = 1 : size(y_ampl, 2)
        for j = 1 : i
            y_ampl(1, i) = y_ampl(1, i) + x(j);
        end
        y_ampl(1, i) = y_ampl(1, i) / i;
    end
    y_per = 0; % TODO

    % set other parameters (mutation probabilities)
    ppar = 0;
    pd = 0;
    pdeg = 0;
    pgen = 0;
    preg = 0;
    po = 0;
    pg = 0;
    plm = 0;
    dl = 0;
    de = 0;

    for iter = 1 : max_iterations
        C = zeros(1, pop_size);
        for i = 1 : pop_size
            mutation(pop_array{i}, ppar, pd, pdeg, pgen, preg, po, pg, plm, dl, de);
            setGlobalx(pop_array{i});
            [~, y] = ode15s(@model_complete, [0, 10 - 0.01], p_con{i});
            % zadnja dva stolpca v y sta resitvi?
            C(1, i) = cost(y(:, end - 1), y_ampl) + cost(y(:, end), y_per);
        end
        % select best
    end
    %select best motif
end