% dP ... diferencial proteinov
function dP = model_complete(t, P)
    G = getGlobalx;
    
    % m ... stevilo genov v enacbi
    m = size(G, 1);
    dP = zeros(1, size(P, 1));

    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%% STRUKTURA GENA %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    % G(i,1)  ... nacin izrazanja
    %             (0 = gensko, 1 = linearna mod., 2 = encimska mod)
    % G(i,2)  ... parameter alfa/beta
    % G(i,3)  ... indeks prvega proteina (P2 pri modifikaciji)
    % G(i,4)  ... mejna koncentracija prvega proteina (P2 pri modifikaciji)
    % G(i,5)  ... indeks drugega proteina
    % G(i,6)  ... mejna koncentracija drugega proteina
    % G(i,7)  ... nacin degradacije
    %             (0 = linearna, 1 = encimska, 2 = aktivna)
    % G(i,8)  ... parameter delta
    % G(i,9)  ... koncentracija, pri kateri je razgradnja polovicna
    % G(i,10) ... indeks drugega proteina

    for i = 1 : m
        if G(i, 1) == 0 % gensko izrazanje
            tf1 = 1;
            tf2 = 1;

            if G(i, 3) > 0 % aktivator
                tf1 = heaviside(P(G(i, 3)) - G(i, 4));
            elseif G(i, 3) < 0 % represor
                tf1 = heaviside(G(i, 4) - P(-G(i, 3)));
            end

            if G(i, 5) > 0 % aktivator
                tf2 = heaviside(P(G(i, 5)) - G(i, 6));
            elseif G(i, 5) < 0 % represor
                tf2 = heaviside(G(i, 6) - P(-G(i, 5)));
            end

            dP(i) = dP(i) + G(i, 2) * tf1 * tf2;
        elseif G(i, 1) == 1 % linearna modifikacija
            % Koncentracija trenutnega proteina se zmanjsa.
            dP(i) = dP(i) - G(i, 2) * P(i);
            dP(G(i, 3)) = dP(G(i, 3)) + G(i, 2) * P(i);
        elseif G(i, 1) == 2 % encimska modifikacija
            dP(i) = dP(i) - G(i, 2) * (P(i) / (G(i, 4) + P(i)));
            dP(G(i, 3)) = dP(G(i, 3)) + G(i, 2) * (P(i) / (G(i, 4) + P(i)));
        end

        if(i ~= 1)
            if G(i, 7) == 0 % linearna degradacija
                dP(i) = dP(i) - G(i, 8) * P(i);
            elseif G(i, 7) == 1 % encimska degradacija
                dP(i) = dP(i) - G(i, 8) * P(i) / (G(i, 9) + P(i));
            elseif G(i, 7) == 2 % aktivna degradacija
                dP(i) = dP(i) - G(i, 8) * P(i) * P(G(i, 10));
            end
        end
    end

    Par = getGlobalAP();
    dP(1) = Par(1) * ((2 * pi) / Par(2)) * cos((2 * t * pi) / Par(2));
    dP = dP';
end
