function NG = mutation(G, pd, po, pg, plm, dl, de)
    % G   ... gensko regulatorno omrezje
    % pd  ... verjetnost dodajanja proteina
    % po  ... verjetnost odstranitve proteina
    % pg  ... verjetnost genskega izrazanja
    % plm ... verjetnost linearne modifikacije
    % dl  ... verjetnost linearne degradacije
    % de  ... verjetnost encimske degradacije
    sP = size(G, 1); % stevilo proteinov

    % TODO: Verjetno ne zelimo, da se tako veliko sprememb zgodi.
    %       Najboljse je, da tudi nakljucno izberemo eno izmed njih.
    % TODO: Verjetno ne zelimo izbrisati izhodnih proteinov. Ali naj bosta
    %       izhodna proteina vnaprej dolocena ali sta vedno zadnja proteina
    %       izhodna?

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sprememba obstojecega parametra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sel = ceil(rand() * (sP - 1) + 1); % ne mutiramo prvega gena
    params = [2, 4, 6, 8, 9];
    idx = params(ceil(rand() * length(params)));
    G(sel, idx) = G(sel, idx) * (rand() * 2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dodajanje novega proteina %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(sP < 10 && rand() <= pd)
        P = zeros(1, 10);

        % TODO: Kako naj se nastavi P(2)? Trenutno dobi vrednost med 1 in 100.
        % TODO: Enako za P(4), P(6), P(8) in P(9) ter drugje.
        P(2) = rand() * 100 + 1;
        P(8) = rand() * 100 + 1;

        % nacin izrazanja
        res = rand();
        if(res < pg) % gensko izrazanje
            P(1) = 0;

            if(rand() >= 0.5) % prvi transkripcijski faktor
                P(3) = randsample(setdiff(1 : sP, sel), 1);
                P(4) = rand() * 100 + 1;
                if(rand() < 0.5)
                    P(3) = -P(3);
                end
            end

            if(rand >= 0.5) % drugi transkripcijski faktor
                P(5) = randsample(setdiff(1 : sP, [abs(P(3)), sel]), 1);
                P(6) = rand() * 100 + 1;
                if(rand() < 0.5)
                    P(5) = -P(5);
                end
            end
        elseif(res < pg + plm) % linearna modifikacija
            P(1) = 1;
            P(3) = randsample(setdiff(1 : sP, sel), 1);
        else % encimska modifikacija
            P(1) = 2;
            P(3) = randsample(setdiff(1 : sP, sel), 1);
            P(4) = rand() * 100 + 1;
        end

        % nacin degradacije
        res = rand();
        if(res < dl) % linearna degradacija
            P(7) = 0;
        elseif(res < dl + de) % encimska degradacija
            P(7) = 1;
            P(9) = rand() * 100 + 1;
        else % aktivna degradacija
            P(7) = 2;
            P(10) = randsample(setdiff(1 : sP, sel), 1);
        end

        G = [G; P];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spreminjanje nacina degradacije %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO: Ne upostevaj novo dodanega proteina, ce je bil dodan.
    sel = ceil(rand() * (sP - 1) + 1);
    G(sel, 8) = rand() * 100 + 1;

    res = rand();
    if(res < dl) % linearna degradacija
        G(sel, 7) = 0;
    elseif(res < dl + de) % encimska degradacija
        G(sel, 7) = 1;
        G(sel, 9) = rand() * 100 + 1;
    else % aktivna degradacija
        G(sel, 7) = 2;
        G(sel, 10) = randsample(setdiff(1 : sP, sel), 1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spreminanje nacina generiranja proteina %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO: Ne upostevaj novo dodanega proteina, ce je bil dodan.
    sel = ceil(rand() * (sP - 1) + 1);
    G(sel, 2) = rand() * 100 + 1;

    res = rand();
    if(res < pg) % gensko izrazanje
        G(sel, 1) = 0;

        if(rand() >= 0.5) % prvi transkripcijski faktor
            G(sel, 3) = randsample(setdiff(1 : sP, sel), 1);
            G(sel, 4) = rand() * 100 + 1;
            if(rand() < 0.5)
                G(sel, 3) = -G(sel, 3);
            end
        end

        if(rand >= 0.5) % drugi transkripcijski faktor
            G(sel, 5) = randsample(setdiff(1 : sP, [abs(G(sel, 3)), sel]), 1);
            G(sel, 6) = rand() * 100 + 1;
            if(rand() < 0.5)
                G(sel, 5) = -G(sel, 5);
            end
        end
    elseif(res < pg + plm) % linearna modifikacija
        G(sel, 1) = 1;
        G(sel, 3) = randsample(setdiff(1 : sP, sel), 1);
    else % encimska modifikacija
        G(sel, 1) = 2;
        G(sel, 3) = randsample(setdiff(1 : sP, sel), 1);
        G(sel, 4) = rand() * 100 + 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spreminjanje regulatorjev %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sel = ceil(rand() * (sP - 1) + 1);
    if(G(sel, 1) == 0)
        % TODO: Nakljucno izberi vrsto spremembe.
        % dodajanje regulatorja
        if(G(sel, 3) == 0)
            G(sel, 3) = randsample(setdiff(1 : sP, [abs(G(sel, 5)), sel]), 1);
        elseif(G(sel, 5) == 0)
            G(sel, 5) = randsample(setdiff(1 : sP, [abs(G(sel, 3)), sel]), 1);
        end

        % brisanje regulatorja
        % TODO: Lahko odstranimo regulator, tudi ce je samo eden
        %       -> brezpogojno izrazanje.
        if(G(sel, 3) ~= 0 && G(sel, 5) ~= 0)
            G(sel, randsample([3, 5], 1)) = 0;
        end

        % spreminjanje vrste regulatorja
        if(rand() >= 0.5)
            G(sel, 3) = -G(sel, 3);
        end
        if(rand() >= 0.5)
            G(sel, 5) = -G(sel, 5);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%
    % odstranitev proteina %
    %%%%%%%%%%%%%%%%%%%%%%%%
    % TODO: Namesto spremembe lahko parametre tudi odstranimo
    %       (pri tem je potrebno spremeniti nacin npr. degradacije).
    if(sP > 3 && rand() < po)
        sel = ceil(rand() * (sP - 1) + 1);

        for i = 1 : sP
            if sel ~= i
                if(abs(G(i, 3)) == sel)
                    G(sel, 3) = randsample(setdiff(1 : sP, [sel, abs(G(i, 5)), i]), 1);
                    if(G(i,1) == 0 && rand() >= 0.5)
                        G(sel, 3) = -G(sel, 3);
                    end
                end

                if(abs(G(i, 5)) == sel)
                    G(sel, 5) = randsample(setdiff(1 : sP, [sel, abs(G(i, 3)), i]), 1);
                    if(rand() >= 0.5)
                        G(sel, 5) = -G(sel,5);
                    end
                end

                if(abs(G(i, 10)) == sel)
                    G(sel, 10) = randsample(setdiff(1 : sP, [sel, i]), 1);
                end
            end
        end

        [G, ~] = removerows(G, 'ind', sel);
    end

    NG = G;
end