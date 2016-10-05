function [winner summary] = play_TexasHoldem(handcards, nplayers, N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulation of many Texas Hold'em Poker games to estimate winning
% probability of given handcards and given number of players.
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

% check input
if nplayers < 2 || nplayers > 23 || int32(nplayers) ~= nplayers
    error('number of players must be an integer (minimum 2, maximum 23).')
end

if length(handcards) ~= 2 || ~all(ismember(handcards,cards)) ...
   || strcmp(handcards{1},handcards{2})
    error('choose two valid handcards.')
end

if N < 1 || N > 1e6
    warning('N too small or too big... select N=1000')
    N = 1000;
end

nwins = 0;
ntie = 0;

winner = struct('index', {}, 'handrank', {}, 'handname', {}, 'hand', {});
    
% wb = waitbar(0,'Playing Poker...');
tic
for i = 1:N
    
    % distribute random handcards to other players
    % (one by one; from remaining stack)
    playercards = cell(nplayers,2);
    commoncards = cell(1,5);
    
    playercards(1,:) = handcards;
    
    remainstack = cards(~ismember(cards,playercards(1,:)));
    
    for c = 1:2
        % distribute hand cards one by one
        for n = 2:nplayers
            icard = randi(length(remainstack),[1 1]);
            playercards(n,c) = remainstack(icard);
            remainstack = remainstack(~ismember(remainstack,playercards(n,c)));
        end
    end
    
    % randomly select five common cards (3 flop, 1 turn and 1 river)
    % (one by one; from remaining stack)
    for c = 1:5
        icard = randi(length(remainstack),[1 1]);
        commoncards(c) = remainstack(icard);
        remainstack = remainstack(~ismember(remainstack,commoncards(c)));
    end
      
    % check hand ranking for each player
    for n = 1:nplayers
        pokerHand(n) = check_PokerHand(playercards(n,:),commoncards);
        handrank(n) = pokerHand(n).rank;
    end
    
    % the player with highest handrank wins
    ihighest = find(handrank == max(handrank));
    
    if length(ihighest) == 1
        % we have a clear winner
        winner(i).index = ihighest;
        winner(i).handrank = handrank(ihighest);
        winner(i).handname = pokerHand(ihighest).name;
        winner(i).hand = pokerHand(ihighest).cards;
        if winner(i).index == 1
            % player has won
            nwins = nwins + 1;
        end
    else
        % multiple players have same rank
        % find winner or split pot
        winmatrix = zeros(length(ihighest),length(ihighest));
        for m = 1:length(ihighest)
            for n = m+1:length(ihighest)
                winmatrix(m,n) = compare_TieHands(pokerHand(ihighest(m)),pokerHand(ihighest(n)));
                if winmatrix(m,n) == 1
                    winmatrix(n,m) = 0;
                elseif winmatrix(m,n) == 2
                    winmatrix(m,n) = 0;
                    winmatrix(n,m) = 1;
                end
            end
        end
        iwinner = find(sum(winmatrix,2) == max(sum(winmatrix,2)));
        if length(iwinner) == 1
            winner(i).index = iwinner;
            winner(i).handrank = handrank(iwinner);
            winner(i).handname = pokerHand(iwinner).name;
            winner(i).hand = pokerHand(iwinner).cards;
            if winner(i).index == 1
                % player has won
                nwins = nwins + 1;
            end
        elseif ~isempty(find(iwinner == 1, 1))
            % split pot (player gets share)
            ntie = ntie + 1;
            winner(i).index = 0;
            winner(i).handrank = handrank(iwinner(1));
            winner(i).handname = pokerHand(iwinner(1)).name;
            winner(i).hand = 'multiple';
        else
            % split pot (player does not get share)
            winner(i).index = 0;
            winner(i).handrank = handrank(iwinner(1));
            winner(i).handname = pokerHand(iwinner(1)).name;
            winner(i).hand = 'multiple';
        end
    end
    
%     waitbar(i/N,wb)
 
end
toc
% close(wb)

% make analyses
summary.prob_win = roundn(nwins*100/N, -1);
summary.prob_tie = roundn(ntie*100/N, -1);

disp([num2str(summary.prob_win) '% wins'])
disp([num2str(summary.prob_tie) '% split pot share'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%