function entropy = entroInf(mean_vec_cat,var_vec_cat)
    entropy = [];
    for i = 1:size(mean_vec_cat,2)
        mean_vec_cat_undertest = mean_vec_cat(:,i);
        var_vec_cat_undertest = var_vec_cat(:,i);
        vector_inf = var_vec_cat_undertest./mean_vec_cat_undertest;
    
        % Generate the probability distribution
        [counts_cat, ~] = histcounts(vector_inf,'Normalization', 'probability');
        
        % Remove zero probabilities (as log(0) is undefined)
        nonzero_probs_cat = counts_cat(counts_cat > 0);

        % Compute entropy
        entropy = [entropy , -sum(nonzero_probs_cat .* log2(nonzero_probs_cat))];
    end
    
end