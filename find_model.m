function find_model
    %number of iterations
    max_iterations = 1000;
    %population size 
    pop_size = 20;
    % create initaial population
    pop_array = cell(1,pop_size); %cell array 
    %creates cell array of empty matrices
    p_con = cell(1,pop_size); 
    %protein concentrations
    for i = 1:pop_size
        %zacetno stevilo proteinov
        M = ones(3,10); 
        M(:,1) = 0; %spremeni tip izrazanja na gensko izrazanje
        M(:,5) = 0;
        M(:,6) = 0;
        M(:,7) = 0; %spremeni tip degradacije na linearno alfa = 1
        pop_array{i} = M;
        
        p_con{i} = zeros(1,3);        
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

    for iter = 1:max_iterations
        % mutate 
        % evaluate on all learning examples (cost function)
        % select best
    end
    %select best motif
end