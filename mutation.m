function NG = mutation(G, ppar, pd, pdeg, pgen, preg, po, pg, plm, dl, de)
    % G    ... gensko regulatorno omrezje
    % pd   ... verjetnost dodajanja proteina
    % po   ... verjetnost odstranitve proteina
    % pg   ... verjetnost genskega izrazanja
    % plm  ... verjetnost linearne modifikacije
    % dl   ... verjetnost linearne degradacije
    % de   ... verjetnost encimske degradacije
    % ppar ... verjetnost spremembe parametrov
    % pdeg ... verjetnost spremembe degradacije
    % pgen ... verjetnost spremembe generiranja proteina
    % preg ... verjetnost spremembe regulatorjev

    sP = size(G, 1); % stevilo proteinov

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sprememba obstojecega parametra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (rand() <= ppar)
        sel = ceil(rand() * (sP - 1) + 1); % ne mutiramo prvega gena
        params = [2, 4, 6, 8, 9];
        idx = params(ceil(rand() * length(params)));
        G(sel, idx) = G(sel, idx) * (rand() * 2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dodajanje novega proteina %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(sP < 10 && rand() <= pd)
        P = zeros(1, 10);

        P(2) = 1;
        P(8) = 1;

        res = rand(); % nacin izrazanja
        if(res < pg) % gensko izrazanje
            P(1) = 0;

            if(rand() > 0.5) % prvi transkripcijski faktor
                P(3) = randsample(1 : (sP + 1), 1);
                P(4) = 1;
                if(rand() < 0.5)
                    P(3) = -P(3);
                end
            end

            if(rand() > 0.5) % drugi transkripcijski faktor
                P(5) = randsample(setdiff(1 : (sP + 1), abs(P(3))), 1);
                P(6) = 1;
                if(rand() < 0.5)
                    P(5) = -P(5);
                end
            end
        elseif(res < pg + plm) % linearna modifikacija
            P(1) = 1;
            P(3) = randsample(1 : sP, 1);
        else % encimska modifikacija
            P(1) = 2;
            P(3) = randsample(1 : sP, 1);
            P(4) = 1;
        end

        res = rand(); % nacin degradacije
        if(res < dl) % linearna degradacija
            P(7) = 0;
        elseif(res < dl + de) % encimska degradacija
            P(7) = 1;
            P(9) = 1;
        else % aktivna degradacija
            P(7) = 2;
            P(10) = randsample(1 : sP, 1);
        end

        G = [G; P];
        sP = sP + 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spreminjanje nacina degradacije %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(rand() <= pdeg)
        sel = ceil(rand() * (sP - 1) + 1);
        G(sel, 8) = 1;

        res = rand();
        if(res < dl) % linearna degradacija
            G(sel, 7) = 0;
        elseif(res < dl + de) % encimska degradacija
            G(sel, 7) = 1;
            G(sel, 9) = 1;
        else % aktivna degradacija
            G(sel, 7) = 2;
            G(sel, 10) = randsample(setdiff(1 : sP, sel), 1);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spreminjanje nacina generiranja proteina %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(rand() <= pgen)
        sel = ceil(rand() * (sP - 1) + 1);
        G(sel, 2) = 1;

        res = rand();
        if(res < pg) % gensko izrazanje
            G(sel, 1) = 0;

            if(rand() > 0.5) % prvi transkripcijski faktor
                G(sel, 3) = randsample(setdiff(1 : sP, abs(G(sel, 5))), 1);
                G(sel, 4) = 1;
                if(rand() < 0.5)
                    G(sel, 3) = -G(sel, 3);
                end
            end

            if(rand() > 0.5) % drugi transkripcijski faktor
                G(sel, 5) = randsample(setdiff(1 : sP, abs(G(sel, 3))), 1);
                G(sel, 6) = 1;
                if(rand() < 0.5)
                    G(sel, 5) = -G(sel, 5);
                end
            end
        elseif(res < pg + plm) % linearna modifikacija
            G(sel, 1) = 1;
            G(sel, 3) = randsample(setdiff(1 : sP, sel), 1);
            G(sel, 5) = 0;
        else % encimska modifikacija
            G(sel, 1) = 2;
            G(sel, 3) = randsample(setdiff(1 : sP, sel), 1);
            G(sel, 4) = 1;
            G(sel, 5) = 0;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spreminjanje regulatorjev %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(rand() <= preg)
        sel = ceil(rand() * (sP - 1) + 1);
        if(G(sel, 1) == 0)
            res = rand();
            if(res <= 0.30 && (G(sel, 3) == 0 || G(sel, 5) == 0)) % dodajanje regulatorja
                if(G(sel, 3) == 0)
                    G(sel, 3) = randsample(setdiff(1 : sP, abs(G(sel, 5))), 1);
                    G(sel, 4) = 1;
                elseif(G(sel, 5) == 0)
                    G(sel, 5) = randsample(setdiff(1 : sP, abs(G(sel, 3))), 1);
                    G(sel, 6) = 1;
                end
            elseif(res > 0.70 && (G(sel, 3) ~= 0 || G(sel, 5) ~= 0)) % brisanje regulatorja
                if(rand() > 0.5 && G(sel, 3) ~= 0)
                    G(sel, 3) = 0;
                else
                    G(sel, 5) = 0;
                end
            elseif(G(sel, 3) ~= 0) % spreminjanje vrste prvega regulatorja
                G(sel, 3) = -G(sel, 3);
            else % spreminjanje vrste drugega regulatorja
                G(sel, 5) = -G(sel, 5);
            end

            if(rand() > 0.5 && G(sel, 3) ~= 0) % sprememba koncentracije
                G(sel, 4) = G(sel, 4) * (rand() * 2);
            elseif(G(sel, 5) ~= 0)
                G(sel, 6) = G(sel, 6) * (rand() * 2);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%
    % odstranitev proteina %
    %%%%%%%%%%%%%%%%%%%%%%%%
    if(sP > 3 && rand() < po)
        sel = ceil(rand() * (sP - 1) + 1);
        G(sel, :) = [];
        sP = sP - 1;

        for i = 1 : sP
            if(G(i, 1) == 0)
                if(abs(G(i, 3)) == sel)
                    G(i, 3) = randsample(setdiff(1 : sP, abs(G(i, 5))), 1);
                    if(rand() > 0.5)
                        G(i, 3) = -G(i, 3);
                    end
                end

                if(abs(G(i, 5)) == sel)
                    G(i, 5) = randsample(setdiff(1 : sP, abs(G(i, 3))), 1);
                    if(rand() > 0.5)
                        G(i, 5) = -G(i, 5);
                    end
                end
            else
                if(abs(G(i, 3)) == sel)
                    G(i, 3) = randsample(setdiff(1 : sP, i), 1);
                end
                G(i, 5) = 0;
            end

            if(abs(G(i, 10)) == sel)
                G(i, 10) = randsample(setdiff(1 : sP, i), 1);
            end

            if(G(i, 3) > sel)
                G(i, 3) = G(i, 3) - 1;
            elseif(abs(G(i, 3)) > sel)
                G(i, 3) = G(i, 3) + 1;
            end

            if(G(i, 5) > sel)
                G(i, 5) = G(i, 5) - 1;
            elseif(abs(G(i, 5)) > sel)
                G(i, 5) = G(i, 5) + 1;
            end

            if(G(i, 10) > sel)
                G(i, 10) = G(i, 10) - 1;
            end
        end
    end

    NG = G;
end
