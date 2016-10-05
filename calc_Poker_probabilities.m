P = 4;

ranks = {'2','3','4','5','6','7','8','9','10','j','q','k','a'};
nranks = length(ranks);

for i = 1:nranks
    disp(i)
    for j = i:nranks
        % combos with different colors
        cards_diffcol{i,j} = {['c' ranks{i}], ['s' ranks{j}]};
        [~, prob_diffcol(i,j)] = play_TexasHoldem(cards_diffcol{i,j}, P, 10000);
        
        % combos with the same color
        if j ~= i
            cards_samecol{i,j} = {['s' ranks{i}], ['s' ranks{j}]};
            [~, prob_samecol(i,j)] = play_TexasHoldem(cards_samecol{i,j}, P, 10000);
        end
    end
end

% save('/home/kay/Matlab/Poker_probabilities/prob_8p_10k.mat','*_8p_10k')

pwin_diffcol_4p_10k = nan(13,13);
ptie_diffcol_4p_10k = nan(13,13);
pwin_samecol_4p_10k = nan(12,13);
ptie_samecol_4p_10k = nan(12,13);
cards_diffcol_4p_10k = cell(13,13);
cards_samecol_4p_10k = cell(12,13);

for i = 1:13
    for j = 1:13
        if ~isempty(prob_diffcol(i,j).prob_win)
            pwin_diffcol_4p_10k(i,j) = prob_diffcol(i,j).prob_win;
            ptie_diffcol_4p_10k(i,j) = prob_diffcol(i,j).prob_tie;
            dum = cards_diffcol{i,j}(2);
            cards_diffcol_4p_10k{i,j} = dum{:}(2:end);
        end
        
        if i < 13 && ~isempty(prob_samecol(i,j).prob_win)
            pwin_samecol_4p_10k(i,j) = prob_samecol(i,j).prob_win;
            ptie_samecol_4p_10k(i,j) = prob_samecol(i,j).prob_tie;
            dum = cards_samecol{i,j}(2);
            cards_samecol_4p_10k{i,j} = dum{:}(2:end);
        end
    end
end

save('/home/kay/Matlab/Poker_probabilities/prob_4p_10k.mat','*_4p_10k')