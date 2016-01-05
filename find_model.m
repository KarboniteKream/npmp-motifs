function find_model
    % number of iterations
    max_iterations = 1000;
    % population size
    pop_size = 20;
    
    % create initial population
    % creates cell array of empty matrices
    pop_array = cell(1, pop_size * 2); %cell array
    % protein concentrations
    p_con = cell(1, pop_size * 2);
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
    % TODO: Dodaj vhodni signal v osebke kot prvi protein.
    t = 0 : 0.01 : 9.99;
    f = 10;
    x = sin(f * t);
    y_ampl = mean(x) * ones(size(x));
    y_per = ((2 * pi) / f) * ones(size(x));

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

    % cost posameznih osebkov
    C = zeros(1, pop_size * 2);

    % glavna zanka
    for iter = 1 : max_iterations
        % mutacija - dodamo 20 novih osebkov
        for i = 1 : pop_size
            pop_array{20 + i} = mutation(pop_array{i}, ppar, pd, pdeg, pgen, preg, po, pg, plm, dl, de);
        end

        % simulacija
        for i = 1 : (pop_size * 2)
            setGlobalx(pop_array{i});
            [~, y] = ode15s(@model_complete, [0, 10 - 0.01], p_con{i});
            C(1, i) = cost(y(:, 2), y_ampl) + cost(y(:, 3), y_per);
        end

        % select best
        [~, sort_idx] = sort(C);
        pop_array{1 : pop_size} = pop_array{sort_idx(1 : pop_size)};
    end

    % select best motif
    pop_array{1}
end