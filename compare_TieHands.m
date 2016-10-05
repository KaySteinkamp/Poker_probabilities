function iwin = compare_TieHands(hand1,hand2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compares two Poker Hands of equal rank (e.g. two Full Houses) and
% determines the winner or split pot.
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
if ~isfield(hand1,'cards') || ~isfield(hand2,'cards') ||...
   ~isfield(hand1,'rank') || ~isfield(hand2,'rank') 
    error('Invalid input structures')
end

% initialize iwin to split pot
iwin = 0;

% compare hands
if ~isequal(hand1.rank,hand2.rank)
    [winrank, iwin] = max([hand1.rank,hand2.rank]);
    disp(['hands have different ranks; clear winner is ' ...
          num2str(iwin) ' with ' hranks(winrank)])
else
    % sort card ranks in descending order
    % not done any more
%     r1 = sort(hand1.cardranks,1,'descend');
%     r2 = sort(hand2.cardranks,1,'descend'); 
    r1 = hand1.cardranks;
    r2 = hand2.cardranks;   
    switch hand1.rank
        
        case 1
            % High Card
            i = find(r1 - r2, 1, 'first');
            if ~isempty(i)
                [~, iwin] = max([r1(i),r2(i)]);
            else
                % Split Pot
                iwin = 0;
            end
              
            
        case 2
            % Pair
            nc1 = zeros(length(ranks),1);
            nc2 = zeros(length(ranks),1);
            for r = 1:length(ranks)
                nc1(r) = length(find(r1 == r));
                nc2(r) = length(find(r2 == r));
            end
            p1 = find(nc1 == 2, 1, 'last');
            p2 = find(nc2 == 2, 1, 'last');
            
            if p1 ~= p2
                [~, iwin] = max([p1,p2]);
            else
                % both have the same Pair
                rr1 = r1(r1 ~= p1);
                rr2 = r2(r2 ~= p2);
                
                i = find(rr1 - rr2, 1, 'first');
                if ~isempty(i)
                    [~, iwin] = max([rr1(i),rr2(i)]);
                else
                    % Split Pot
                    iwin = 0; 
                end
            end
            
            
        case 3
            % Two Pairs
            nc1 = zeros(length(ranks),1);
            nc2 = zeros(length(ranks),1);
            for r = 1:length(ranks)
                nc1(r) = length(find(r1 == r));
                nc2(r) = length(find(r2 == r));
            end
            p11 = find(nc1 == 2, 1, 'first');
            p12 = find(nc1 == 2, 1, 'last');
            p21 = find(nc2 == 2, 1, 'first');
            p22 = find(nc2 == 2, 1, 'last');
            
            if p12 ~= p22
                % high pairs different
                [~, iwin] = max([p12,p22]);
            elseif p11 ~= p21
                % small pairs different
                [~, iwin] = max([p11,p21]);
            else
                % both pairs equal
                % check for remaining high card
                rr1 = r1(r1 ~= p11 & r1 ~= p12);
                rr2 = r2(r2 ~= p21 & r2 ~= p22);
                if ~isequal(rr1,rr2)
                    [~, iwin] = max([rr1,rr2]);
                else
                    % Split Pot
                    iwin = 0; 
                end
            end
                
                
        case 4
            % Three of a kind
            nc1 = zeros(length(ranks),1);
            nc2 = zeros(length(ranks),1);
            for r = 1:length(ranks)
                nc1(r) = length(find(r1 == r));
                nc2(r) = length(find(r2 == r));
            end
            t1 = find(nc1 == 3, 1, 'last');
            t2 = find(nc2 == 3, 1, 'last');
            
            if t1 ~= t2
                [~, iwin] = max([t1,t2]);
            else
                % both have the same Three
                rr1 = r1(r1 ~= t1);
                rr2 = r2(r2 ~= t2);
                
                i = find(rr1 - rr2, 1, 'first');
                if ~isempty(i)
                    [~, iwin] = max([rr1(i),rr2(i)]);
                else
                    % Split Pot
                    iwin = 0; 
                end
            end
            
            
        case 5
            % Straight
            % Note that Straights are sorted in ascending order
            if isequal(r1(end),r2(end))
                % Split Pot
                iwin = 0;
            else
                [~, iwin] = max([r1(end),r2(end)]);
            end
            
            
        case 6
            % Flush
            % Note that Flushes are sorted in ascending order
            i = find(r1 - r2, 1, 'last');
            if ~isempty(i)
                [~, iwin] = max([r1(i),r2(i)]);
            else
                % Split Pot
                iwin = 0;
            end
            
            
        case 7
            % Full House
            nc1 = zeros(length(ranks),1);
            nc2 = zeros(length(ranks),1);
            for r = 1:length(ranks)
                nc1(r) = length(find(r1 == r));
                nc2(r) = length(find(r2 == r));
            end
            t1 = find(nc1 == 3, 1, 'last');
            p1 = find(nc1 == 2, 1, 'last');
            t2 = find(nc2 == 3, 1, 'last');
            p2 = find(nc2 == 2, 1, 'last');
            
            if t1 ~= t2
                [~, iwin] = max([t1,t2]);
            elseif p1 ~= p2
                % both have the same Three
                % check for higher Pair
                [~, iwin] = max([p1,p2]);
            else
                % both have the same Three and same Pair,
                % so Split Pot
                iwin = 0;
            end
                            
            
        case 8
            % Four of a kind
            nc1 = zeros(length(ranks),1);
            nc2 = zeros(length(ranks),1);
            for r = 1:length(ranks)
                nc1(r) = length(find(r1 == r));
                nc2(r) = length(find(r2 == r));
            end
            f1 = find(nc1 == 4, 1, 'last');
            f2 = find(nc2 == 4, 1, 'last');
            
            if f1 ~= f2
                [~, iwin] = max([f1,f2]);
            else
                % both have the same Four
                % check for remaining high card
                rr1 = r1(r1 ~= f1);
                rr2 = r2(r2 ~= f2);
                if ~isequal(rr1,rr2)
                    [~, iwin] = max([rr1,rr2]);
                else
                    % Split Pot
                    iwin = 0; 
                end
            end
            
            
        case 9
            % Straight Flush   
            % Note that Straights are sorted in ascending order
            if isequal(r1(end),r2(end))
                % Split Pot
                iwin = 0;
            else
                [~, iwin] = max([r1(end),r2(end)]);
            end
            
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%