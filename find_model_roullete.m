function find_model_roullete
    % number of iterations
    max_iterations = 500;
    % population size
    pop_size = 10;
    % stevilo ucnih signalov
    signal_number = 5;

    % create initial population
    % creates cell array of empty matrices
    pop_array = cell(1, pop_size * 2); %cell array
    % protein concentrations
    for i = 1 : pop_size
        % zacetno stevilo proteinov
        M = ones(3, 10);
        M(:, 1) = 0; % spremeni tip izrazanja na gensko izrazanje
        M(:, 5) = 0;
        M(:, 6) = 0;
        M(:, 7) = 0; % spremeni tip degradacije na linearno alfa = 1
        pop_array{i} = M;
    end

    % create learning examples
    S = zeros(signal_number,2);
    for i = 1:signal_number
        S(i,1) = randsample(1 : 5, 1); % predstavlja amplitudo
        S(i,2) = randsample(1 : 4, 1); % predstavlja perido
    end
    % TODO: Dodaj vhodni signal v osebke kot prvi protein.

    % set other parameters (mutation probabilities)
    ppar = 0.85;   % ppar ... verjetnost spremembe parametrov
    pd   = 0.40;   % pd   ... verjetnost dodajanja proteina
    pdeg = 0.35;   % pdeg ... verjetnost spremembe degradacije  
    pgen = 0.25; % pgen ... verjetnost spremembe generiranja proteina
    preg = 0.40;  % preg ... verjetnost spremembe regulatorjev
    po   = 0.35; % po   ... verjetnost odstranitve proteina
    pg   = 0.50; % pg   ... verjetnost genskega izrazanja
    plm  = 0.25; % plm  ... verjetnost linearne modifikacije
    dl   = 0.40;  % dl   ... verjetnost linearne degradacije
    de   = 0.35;  % de   ... verjetnost encimske degradacije

    % t = 0 : 0.01 : 9.99;
    t = [0, 20];

    % glavna zanka
    for iter = 1 : max_iterations
        fprintf('iter=%02d', iter);
        % cost posameznih osebkov se ponastavi pri vsaki iteraciji
        C = zeros(1, pop_size * 2);
        % mutacija - dodamo 20 novih osebkov
        for i = 1 : pop_size
            pop_array{pop_size + i} = mutation(pop_array{i}, ppar, pd, pdeg, pgen, preg, po, pg, plm, dl, de);
        end

        % simulacija
        for i = 1 : (pop_size * 2)
            setGlobalx(pop_array{i});
            conc = zeros(1, size(pop_array{i}, 1));
            for j = 1:signal_number
                % nastavi periodo in amplitudo
                setGlobalAP(S(j, :)); % AP as amplitude and period
                % initial protein values are zeros by default
                conc(1) = S(j,1);
                [~, y] = ode45(@model_complete, t, conc);
                C(1, i) = C(1, i) + cost(y(:, 3), S(j, 2)); % + cost(y(:, 2), S(j, 1));
            end
        end

        % select best
        %[~, sort_idx] = sort(C);
        %for i = 1 : pop_size
        %    pop_array{i} = pop_array{sort_idx(i)};
        %end
        fprintf(', cost=%d\n', round(min(C)));
        if round(min(C)) < 50
            break; %terminate search when good solution is found
        end
        %roulete selection
        %total sum of cost
        P = sum(C); %total cost
        F = C;
        F = F/P; %noramlize to get probability
        F = 1./F; %substract to get reverse, biggest cost smaller probability
        P = sum(F); %normalize factor
        F = F/P; %C is now selecct probability
        F = cumsum([0,F]); %comulative sum  
        i = 1; %elitizem always keep first one
        selected = cell(1, pop_size); %cell array
        while i < pop_size 
            i = i + 1;
            sel = rand(); %selector from 0 to 1
            j = sum(sel >= F);    
            selected{i} = pop_array{j};
        end
        for i = 2 : pop_size
            pop_array{i} = selected{i}; %revrite to first 10 of population
        end
    end

    % select best motif
    [~, sort_idx] = sort(C); %C is now a probability
    best = pop_array{sort_idx(1)}
    fprintf('\namplituda=%d, perioda=%d\n', S(1, 1), S(1, 2));
    setGlobalx(best);
    conc = zeros(1, size(best, 1));
    for i = 1:size(S,1)
        figure(i);
        setGlobalAP(S(i, :));
        conc(1) = S(i,1);
        [T, y] = ode45(@model_complete, t, conc);
        plot(T, y(:,[1,3]));
        legend('sinus', 'perioda', 'Location', 'northwest');
    end
    S
end
