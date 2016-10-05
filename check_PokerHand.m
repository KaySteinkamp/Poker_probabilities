function pokerHand = check_PokerHand(handcards, commoncards)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determination of the best Poker hand possible with given handcards (2)
% and given commoncards (5).
% 
% Author: Kay Steinkamp
% Date: Feb 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define color and rank abbreviations
% Kreuz = C(lubs)
% Pik = S(pades)
% Herz = H(earts)
% Karo = D(iamonds)
% ---
% Ace:   A
% King:  K
% Queen: Q
% Jack:  J
% 10..2
colors = {'c','s','h','d'};
ranks = {'2','3','4','5','6','7','8','9','10','j','q','k','a'};

% define valid cards (i.e. combos of color and rank)
ncards = length(colors)*length(ranks);
cards = cell(ncards,1);
for c = 1:length(colors)
    for r = 1:length(ranks)
        i = (c-1)*length(ranks) + r;
        cards{i} = [colors{c} ranks{r}];
    end
end

% define valid hands
hranks = {'High Card','Pair','Two Pairs','Three of a kind','Straight',...
          'Flush','Full House','Four of a kind','Straight Flush'};

% check input
if length(handcards) ~= 2 || ~all(ismember(handcards,cards))
    error('invalid handcards.')
end
if length(commoncards) ~= 5 || ~all(ismember(commoncards,cards))
    error('invalid commoncards.')
end

% all cards must be unique
sevencards = cat(1,reshape(handcards,2,1),reshape(commoncards,5,1));
if length(unique(sevencards)) ~= 7
    error('multiple card(s)!')
end

% decompose sevencards into sevencolors and sevenranks/n
sevencards = char(sevencards);
sevencolors = sevencards(:,1);
sevenranks = sevencards(:,2:end);
sevenrankn = nan(length(sevenranks),1);
for c = 1:7
    sevenrankn(c) = find(strcmp(strtrim(sevenranks(c,:)),ranks));
end
[sevenrankn isort] = sort(sevenrankn);
sevencolors = sevencolors(isort);
sevencards = sevencards(isort,:);

rankdiff = diff(sevenrankn);

pokerHand.rank = 0;

% Identify best hand
% -------------------

% seek for multiple card ranks
if length(sevenrankn) == length(unique(sevenrankn))
    
    % neither Pair nor multiple of a kind
    % there is always a High Card: find it first
    % then check for Straight/Flush/Straight Flush 
    
    ohigh1 = (sevenrankn == max(sevenrankn));
    ihigh1 = find(ohigh1);
    ohigh2 = (sevenrankn == max(sevenrankn(~ohigh1)));
    ihigh2 = find(ohigh2);
    ohigh3 = (sevenrankn == max(sevenrankn(~ohigh1 & ~ohigh2)));
    ihigh3 = find(ohigh3);
    ohigh4 = (sevenrankn == max(sevenrankn(~ohigh1 & ~ohigh2 & ~ohigh3)));
    ihigh4 = find(ohigh4);
    ohigh5 = (sevenrankn == max(sevenrankn(~ohigh1 & ~ohigh2 & ~ohigh3 & ~ohigh4)));
    ihigh5 = find(ohigh5);
    
    pokerHand.rank = 1;
    ihand = [ihigh1; ihigh2; ihigh3; ihigh4; ihigh5];
    
    % check for Straight/Flush/Straight Flush    
    consec = find(rankdiff == 1);
    
    % check consecutive cards, be careful with ace
    if (length(consec) == 4 && (consec(4)-consec(1)) == 3) ||...
       (length(consec) == 5 && (consec(5)-consec(1)) == 4) ||...
       (length(consec) == 6 && (consec(6)-consec(1)) == 5)
        % Straight!
        pokerHand.rank = 5;
        ihand = [consec(end-3:end); consec(end)+1];
        % check for Straight Flush
        if length(unique(sevencolors(ihand))) == 1
            % Straight Flush!
            pokerHand.rank = 9;
        elseif length(consec) == 5
            % smaller straight flush still possible
            dum = [consec(1:4); consec(4)+1];
            if length(unique(sevencolors(dum))) == 1
                % Straight Flush!
                pokerHand.rank = 9;
                ihand = dum;
            end
        elseif length(consec) == 6
            % smaller straight flush still possible
            dum = [consec(1:4); consec(4)+1];
            if length(unique(sevencolors(dum))) == 1
                % Straight Flush!
                pokerHand.rank = 9;
                ihand = dum;
            end
            dum = [consec(2:5); consec(5)+1];
            if length(unique(sevencolors(dum))) == 1
                % Straight Flush!
                pokerHand.rank = 9;
                ihand = dum;
            end
        end
        
    elseif length(consec) >= 3 && sevenrankn(end) == 13
        % maybe still a Straight beginning with ace
        iace = length(sevenrankn);
        if isequal(sevenrankn(1:4), [1;2;3;4])
            % Straight!
            % (beginning with ace)
            pokerHand.rank = 5;
            ihand = [iace; 1; 2; 3; 4];
            % check for Straight Flush
            if length(unique(sevencolors(ihand))) == 1
                % Straight Flush!
                pokerHand.rank = 9;
            end
        end
        
    else
        % check for normal Flush
        for c = 1:length(colors)
            icol = find(strcmp(colors{c},cellstr(sevencolors)));
            if length(icol) >= 5
                % Flush!
                pokerHand.rank = 6;
                % select 5 highest cards for Flush
                % note that sevencards are already sorted
                ihand = icol(end-4:end);                
            end
        end
        
    end

else
    
    % we have multiple cards
    % so check for: 
    % Pair/Two Pairs/Three of a kind/Full House/Four of a kind
    
    counts = zeros(length(ranks),1);
    icount = cell(length(ranks),1);
    for r = 1:length(ranks)
        icount{r} = (sevenrankn == r);
        counts(r) = length(find(icount{r}));
    end
    
    ofour = icount(counts == 4);
    if ~isempty(ofour)
        % Four of a kind!
        % Four of a kind excludes Straight Flush: no need to check
        pokerHand.rank = 8;
        if length(ofour) > 1
            error(['Invalid detection: ' num2str(length(ofour)) ' Fours'])
        end
        % select highest card not part of Four
        ihigh = find(sevenrankn == max(sevenrankn(~ofour{1})));
        ihand = [find(ofour{1}); ihigh(1)];
    else
        nthree = length(find(counts == 3));
        othree = icount(counts == 3);
        if ~isempty(othree)
            % select highest Three (two Threes are possible)
            if nthree == 1
                % single Three:
                % check for Full House
                opair = icount(counts == 2);
                if ~isempty(opair)
                    % Full House!
                    % Full House excludes Straight Flush: no need to check
                    pokerHand.rank = 7;
                    % select highest pair (two pairs are possible)
                    % (simply drop smaller pair)
                    opair = opair(end);
                    ihand = [find(othree{1}); find(opair{1})];
                else
                    % Three of a kind!
                    pokerHand.rank = 4;
                    % select highest 2 cards not part of Three
                    ohigh = (sevenrankn == max(sevenrankn(~othree{1})));
                    ihigh = find(ohigh);
                    ohigh2 = (sevenrankn == max(sevenrankn(~othree{1} & ~ohigh)));
                    ihigh2 = find(ohigh2);
                    ihand = [find(othree{1}); ihigh; ihigh2];
                    
                    % need to check for Straight/Flush/Straight Flush
                    % -----------------------------------------------
                    consec = find(rankdiff == 1);                   
                    % check consecutive cards, be careful with ace
                    % and keep in mind we already have a Three
                    if (length(consec) == 4 && ((consec(4)-consec(1)) == 3 ...
                                                || (consec(4)-consec(1)) == 3+2))
                        % Straight!
                        pokerHand.rank = 5;
                        ihand = [consec(end-3:end); consec(end)+1];
                        % check for Straight Flush
                        ithree = find(othree{1});
                        [iact,~,jh] = intersect(ithree,ihand);
                        iother2 = ithree(ithree ~= iact);
                        % check for each member of Three
                        if length(unique(sevencolors(ihand))) == 1
                            % Straight Flush!
                            pokerHand.rank = 9;
                        else
                            dum = ihand;
                            dum(jh) = iother2(1);
                            if length(unique(sevencolors(dum))) == 1
                                % Straight Flush!
                                pokerHand.rank = 9;
                                ihand = dum;
                            else
                                dum = ihand;
                                dum(jh) = iother2(2);
                                if length(unique(sevencolors(dum))) == 1
                                    % Straight Flush!
                                    pokerHand.rank = 9;
                                    ihand = dum;
                                end
                            end
                        end
                    elseif length(consec) == 3 && sevenrankn(end) == 13
                        % maybe still a Straight beginning with ace
                        iace = length(sevenrankn);
                        i1to4 = [consec(end-2:end); consec(end)+1];
                        if isequal(sevenrankn(i1to4), [1;2;3;4])
                            % Straight!
                            % (beginning with ace)
                            pokerHand.rank = 5;
                            ihand = [iace; i1to4];
                            % check for Straight Flush
                            ithree = find(othree{1});
                            [iact,~,jh] = intersect(ithree,ihand);
                            iother2 = ithree(ithree ~= iact);
                            % check for each member of Three
                            if length(unique(sevencolors(ihand))) == 1
                                % Straight Flush!
                                pokerHand.rank = 9;
                            else 
                                dum = ihand;
                                dum(jh) = iother2(1);
                                if length(unique(sevencolors(dum))) == 1
                                    % Straight Flush!
                                    pokerHand.rank = 9;
                                    ihand = dum;
                                else
                                    dum = ihand;
                                    dum(jh) = iother2(2);
                                    if length(unique(sevencolors(dum))) == 1
                                        % Straight Flush!
                                        pokerHand.rank = 9;
                                        ihand = dum;
                                    end
                                end
                            end 
                        end                       
                    else
                        % check for normal Flush
                        for c = 1:length(colors)
                            icol = find(strcmp(colors{c},cellstr(sevencolors)));
                            if length(icol) >= 5
                                % Flush!
                                pokerHand.rank = 6;
                                % select 5 highest cards for Flush
                                % note that sevencards are already sorted
                                ihand = icol(end-4:end);
                            end
                        end
                        
                    end
                    % -----------------------------------------------
                end
            elseif nthree == 2
                % two Threes:
                % Full House!
                % Full House excludes Straight Flush: no need to check
                % consider lower Three as Pair
                pokerHand.rank = 7;
                ilowthree = find(othree{1});
                ihand = [find(othree{2}); ilowthree(1:2)];
            else
                error(['Invalid detection: ' num2str(nthree) ' Threes'])
            end
            
        else
            % Pair or Two Pairs
            npair = length(find(counts == 2));
            opair = icount(counts == 2);
            
            if npair == 3
                % Two Pairs!
                pokerHand.rank = 3;
                % select the two highest pairs
                % (simply drop smallest pair)
                opair = opair(end-1:end);
                % select highest card not part of Two Pairs
                ihigh = find(sevenrankn == max(sevenrankn(~opair{1} & ~opair{2})));
                ihand = [find(opair{2}); find(opair{1}); ihigh(1)];
                
                % having three pairs excludes Straight(Flush): no need to check
                % but need to check for normal Flush
                % -----------------------------------------------
                for c = 1:length(colors)
                    icol = find(strcmp(colors{c},cellstr(sevencolors)));
                    if length(icol) >= 5
                        % Flush!
                        pokerHand.rank = 6;
                        % select 5 highest cards for Flush
                        % note that sevencards are already sorted
                        ihand = icol(end-4:end);
                    end
                end
                % -----------------------------------------------                
                
            elseif npair == 2
                % Two Pairs!
                pokerHand.rank = 3;
                % select highest card not part of Two Pairs
                ihigh = find(sevenrankn == max(sevenrankn(~opair{1} & ~opair{2})));
                ihand = [find(opair{2}); find(opair{1}); ihigh];
                
                % need to check for Straight/Flush/Straight Flush
                % -----------------------------------------------
                consec = find(rankdiff == 1);                   
                % check consecutive cards, be careful with ace
                % and keep in mind we already have two pairs
                if (length(consec) == 4 && ((consec(4)-consec(1)) == 3 ...
                                            || (consec(4)-consec(1)) == 3+2))
                    % Straight!
                    pokerHand.rank = 5;
                    ihand = [consec(end-3:end); consec(end)+1];
                    % check for Straight Flush
                    ipair1 = find(opair{1});
                    ipair2 = find(opair{2});
                    [iact1,~,jh1] = intersect(ipair1,ihand);
                    [iact2,~,jh2] = intersect(ipair2,ihand);
                    isecond1 = ipair1(ipair1 ~= iact1);
                    isecond2 = ipair2(ipair2 ~= iact2);
                    % check for each member of each Pair
                    if length(unique(sevencolors(ihand))) == 1
                        % Straight Flush!
                        pokerHand.rank = 9;
                    else
                        dum = ihand;
                        dum(jh1) = isecond1;
                        if length(unique(sevencolors(dum))) == 1
                            % Straight Flush!
                            pokerHand.rank = 9;
                            ihand = dum;
                        else
                            dum = ihand;
                            dum(jh2) = isecond2;
                            if length(unique(sevencolors(dum))) == 1
                                % Straight Flush!
                                pokerHand.rank = 9;
                                ihand = dum;
                            else
                                dum = ihand;
                                dum(jh1) = isecond1;
                                dum(jh2) = isecond2;
                                if length(unique(sevencolors(dum))) == 1
                                    % Straight Flush!
                                    pokerHand.rank = 9;
                                    ihand = dum;
                                end
                            end
                        end
                    end
                elseif length(consec) == 3 && sevenrankn(end) == 13
                    % maybe still a Straight beginning with ace
                    iace = length(sevenrankn);
                    i1to4 = [consec(end-2:end); consec(end)+1];
                    if isequal(sevenrankn(i1to4), [1;2;3;4])
                        % Straight!
                        % (beginning with ace)
                        pokerHand.rank = 5;
                        ihand = [iace; i1to4];
                        % check for Straight Flush
                        ipair1 = find(opair{1});
                        ipair2 = find(opair{2});
                        [iact1,~,jh1] = intersect(ipair1,ihand);
                        [iact2,~,jh2] = intersect(ipair2,ihand);
                        isecond1 = ipair1(ipair1 ~= iact1);
                        isecond2 = ipair2(ipair2 ~= iact2);
                        % check for each member of each Pair
                        if length(unique(sevencolors(ihand))) == 1
                            % Straight Flush!
                            pokerHand.rank = 9;
                        else
                            dum = ihand;
                            dum(jh1) = isecond1;
                            if length(unique(sevencolors(dum))) == 1
                                % Straight Flush!
                                pokerHand.rank = 9;
                                ihand = dum;
                            else
                                dum = ihand;
                                dum(jh2) = isecond2;
                                if length(unique(sevencolors(dum))) == 1
                                    % Straight Flush!
                                    pokerHand.rank = 9;
                                    ihand = dum;
                                else
                                    dum = ihand;
                                    dum(jh1) = isecond1;
                                    dum(jh2) = isecond2;
                                    if length(unique(sevencolors(dum))) == 1
                                        % Straight Flush!
                                        pokerHand.rank = 9;
                                        ihand = dum;
                                    end
                                end
                            end
                        end
                    end
                else
                    % check for normal Flush
                    for c = 1:length(colors)
                        icol = find(strcmp(colors{c},cellstr(sevencolors)));
                        if length(icol) >= 5
                            % Flush!
                            pokerHand.rank = 6;
                            % select 5 highest cards for Flush
                            % note that sevencards are already sorted
                            ihand = icol(end-4:end);
                        end
                    end
                    
                end
                % -----------------------------------------------
                          
            elseif npair == 1
                % Pair!
                pokerHand.rank = 2;
                % select highest 3 cards not part of Pair
                ohigh = (sevenrankn == max(sevenrankn(~opair{1})));
                ihigh = find(ohigh);
                ohigh2 = (sevenrankn == max(sevenrankn(~opair{1} & ~ohigh)));
                ihigh2 = find(ohigh2);
                ohigh3 = (sevenrankn == max(sevenrankn(~opair{1} & ~ohigh & ~ohigh2)));
                ihigh3 = find(ohigh3);
                ihand = [find(opair{1}); ihigh; ihigh2; ihigh3];
                
                % need to check for Straight/Flush/Straight Flush
                % -----------------------------------------------
                consec = find(rankdiff == 1);
                % check consecutive cards, be careful with ace
                % and keep in mind we already have one pair
                if (length(consec) == 4 && ((consec(4)-consec(1)) == 3 ...
                                           || (consec(4)-consec(1)) == 3+1))                       
                    % Straight!
                    pokerHand.rank = 5;
                    ihand = [consec; consec(end)+1];
                    
                    % check for Straight Flush
                    % possibilities:
                    % (1) 1 straight, does not include member of pair
                    % (2) 1 straight, includes member of pair
                    ipair = find(opair{1});
                    [iact,~,jh] = intersect(ipair,ihand);
                    
                    % check each possibility
                    if isempty(iact)
                        % possibility (1)
                        if length(unique(sevencolors(ihand))) == 1
                            % Straight Flush!
                            pokerHand.rank = 9;
                        end
                    else
                        % possibility (2)
                        isecond = ipair(ipair ~= iact);
                        if length(unique(sevencolors(ihand))) == 1
                            % Straight Flush!
                            pokerHand.rank = 9;
                        else
                            dum = ihand;
                            dum(jh) = isecond;
                            if length(unique(sevencolors(dum))) == 1
                                % Straight Flush!
                                pokerHand.rank = 9;
                                ihand = dum;
                            end
                        end   
                    end
                
                elseif (length(consec) == 5 && ((consec(5)-consec(1)) == 4 ...
                                               || (consec(5)-consec(1)) == 4+1))                        
                    % Straight!
                    pokerHand.rank = 5;
                    ihand = [consec(end-3:end); consec(end)+1];
                    
                    % check for Straight Flush
                    % possibilities:
                    % (3) 2 straights, one includes member of pair 
                    % (4) 2 straights, both include member of pair
                    ihighS = ihand;
                    ilowS = [consec(1:4); consec(4)+1];
                    ipair = find(opair{1});
                    
                    [iactH,~,jH] = intersect(ipair,ihighS);
                    [iactL,~,jL] = intersect(ipair,ilowS);
                    
                    if isempty(iactH) && isempty(iactL)
                        error('Invalid detection: bad combo of Pair and Straight')
                    end
                                      
                    % check each possibility
                    if isempty(iactH)
                        % possibility (3)
                        if length(unique(sevencolors(ihighS))) == 1
                            % Straight Flush!
                            pokerHand.rank = 9;
                        elseif length(unique(sevencolors(ilowS))) == 1
                            % Straight Flush!
                            pokerHand.rank = 9;
                            ihand = ilowS;
                        else
                            isecond = ipair(ipair ~= iactL);
                            dum = ilowS;
                            dum(jL) = isecond;
                            if length(unique(sevencolors(dum))) == 1
                                % Straight Flush!
                                pokerHand.rank = 9;
                                ihand = dum;
                            end
                        end
                        
                    elseif isempty(iactL)
                        % possibility (3)
                        if length(unique(sevencolors(ihighS))) == 1
                            % Straight Flush!
                            pokerHand.rank = 9;
                        else
                            isecond = ipair(ipair ~= iactH);
                            dum = ihighS;
                            dum(jH) = isecond;
                            if length(unique(sevencolors(dum))) == 1
                                % Straight Flush!
                                pokerHand.rank = 9;
                                ihand = dum;
                            elseif length(unique(sevencolors(ilowS))) == 1
                                % Straight Flush!
                                pokerHand.rank = 9;
                                ihand = ilowS;
                            end
                        end
                        
                    else
                        % possibility (4)
                        isecondH = ipair(ipair ~= iactH);
                        isecondL = ipair(ipair ~= iactL);
                        
                        if length(unique(sevencolors(ihighS))) == 1
                            % Straight Flush!
                            pokerHand.rank = 9;
                        else
                            dum = ihighS;
                            dum(jH) = isecondH;
                            if length(unique(sevencolors(dum))) == 1
                                % Straight Flush!
                                pokerHand.rank = 9;
                                ihand = dum;
                            elseif length(unique(sevencolors(ilowS))) == 1
                                % Straight Flush!
                                pokerHand.rank = 9;
                                ihand = ilowS;
                            else
                                dum = ilowS;
                                dum(jL) = isecondL;
                                if length(unique(sevencolors(dum))) == 1
                                    % Straight Flush!
                                    pokerHand.rank = 9;
                                    ihand = dum;
                                end
                            end
                        end   
                    end
                    
                elseif length(consec) == 3 && sevenrankn(end) == 13
                    % maybe still a Straight beginning with ace
                    iace = length(sevenrankn);
                    i1to4 = [consec(end-2:end); consec(end)+1];
                    if isequal(sevenrankn(i1to4), [1;2;3;4])
                        % Straight!
                        % (beginning with ace)
                        pokerHand.rank = 5;
                        ihand = [iace; i1to4];
                        
                        % check for Straight Flush
                        % possibilities:
                        % (1) straight from ace, does not include member of pair
                        % (2) straight from ace, includes member of pair
                        ipair = find(opair{1});
                        [iact,~,jh] = intersect(ipair,ihand);
                        
                        % check each possibility
                        if isempty(iact)
                            % possibility (1)
                            if length(unique(sevencolors(ihand))) == 1
                                % Straight Flush!
                                pokerHand.rank = 9;
                            end
                        else
                            % possibility (2)
                            isecond = ipair(ipair ~= iact);
                            if length(unique(sevencolors(ihand))) == 1
                                % Straight Flush!
                                pokerHand.rank = 9;
                            else
                                dum = ihand;
                                dum(jh) = isecond;
                                if length(unique(sevencolors(dum))) == 1
                                    % Straight Flush!
                                    pokerHand.rank = 9;
                                    ihand = dum;
                                end
                            end
                        end
                    end
                else
                    % check for normal Flush
                    for c = 1:length(colors)
                        icol = find(strcmp(colors{c},cellstr(sevencolors)));
                        if length(icol) >= 5
                            % Flush!
                            pokerHand.rank = 6;
                            % select 5 highest cards for Flush
                            % note that sevencards are already sorted
                            ihand = icol(end-4:end);
                        end
                    end
                    
                end
                % -----------------------------------------------
                
            else
                error('Failed to detect Pair(s)')
            end
    
        end
        
    end
     
end
   
% Apply detected indices
if ~exist('ihand','var') || length(ihand) ~= 5 || pokerHand.rank == 0
    error('Failed to detect Poker hand.')
end
pokerHand.name = hranks(pokerHand.rank);
pokerHand.cards = sevencards(ihand,:);
pokerHand.cardranks = sevenrankn(ihand);
pokerHand.colors = sevencolors(ihand);
pokerHand.indices = ihand;
pokerHand.set = sevencards;
      
% Test output
% if isequal(reshape(handcards,2,1), {'hj';'hq'}) ||...
%    isequal(reshape(handcards,2,1), {'hq';'hj'})
%     disp(pokerHand.set')
%     disp(pokerHand.cards')
%     disp({'Love'})
% else
%     disp(pokerHand.set')
%     disp(pokerHand.cards')
%     disp(pokerHand.name)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%