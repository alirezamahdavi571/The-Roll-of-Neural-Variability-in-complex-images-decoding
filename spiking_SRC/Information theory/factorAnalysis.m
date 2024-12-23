function factorAnalysis(artifact_Mean_PSTH_data, ...
                             body_Mean_PSTH_data, ...
                             natural_Mean_PSTH_data, ...
                             Face_Mean_PSTH_data, ...
                             var_vec)
% factorAnalysisPart2
% This function combines firing rate data and variance data, performs
% factor analysis, calculates contributions, conducts t-tests, visualizes
% loadings, computes correlations, partial correlations, tries additional
% factors, and performs PCA and nonlinear regression.

    % ---------------------
    %  1) PREPARE THE DATA
    % ---------------------

    % Combine firing rate data from different conditions
    firing_rate_matrix = cat(2, ...
        artifact_Mean_PSTH_data, ...
        body_Mean_PSTH_data, ...
        natural_Mean_PSTH_data, ...
        Face_Mean_PSTH_data);

    % Variance data
    var_matrix = var_vec;

    % Reshape firing rate data:
    %   Rows become (condition * ? ), Columns become time slices
    firing_rate_reshaped = reshape( ...
        firing_rate_matrix, ...
        size(firing_rate_matrix,1) * size(firing_rate_matrix,2), ...
        size(firing_rate_matrix,3));

    % Reshape variance data similarly
    var_reshaped = reshape( ...
        var_matrix, ...
        size(var_matrix,1) * size(var_matrix,2), ...
        size(var_matrix,3));

    % Combine and standardize the data
    combined_data = zscore([firing_rate_reshaped, var_reshaped]);

    % Number of time slices (columns) for the firing_rate_reshaped
    n_time_slices = size(firing_rate_reshaped, 2);

    % --------------------------------
    %  2) PERFORM FACTOR ANALYSIS (2F)
    % --------------------------------

    % Factor Analysis with 2 factors
    [loadings, psi, stats] = factoran(combined_data, 2); % 2 factors

    % Separate loadings for firing rate and variance
    firing_rate_loadings = loadings(1:n_time_slices, :);
    variance_loadings = loadings(n_time_slices+1:end, :);

    % ----------------------
    %  3) VISUALIZE LOADINGS
    % ----------------------
    figure('Name','FA Loadings (2 Factors)','Color','w');
    
    % Firing Rate Loadings
    subplot(1, 2, 1);
    bar(firing_rate_loadings);
    xlabel('Time Slices');
    ylabel('Loadings');
    title('Firing Rate Loadings');

    % Variance Loadings
    subplot(1, 2, 2);
    bar(variance_loadings);
    xlabel('Time Slices');
    ylabel('Loadings');
    title('Variance Loadings');

    % ---------------------------------------
    %  4) CALCULATE CONTRIBUTIONS & DOMINANCE
    % ---------------------------------------
    
    % Sum of absolute loadings across each factor for firing rate
    firing_rate_contribution = sum(abs(firing_rate_loadings), 1);
    % Sum of absolute loadings across each factor for variance
    variance_contribution = sum(abs(variance_loadings), 1);

    % Display contributions
    disp('Firing Rate Contributions (Factor 1, Factor 2):');
    disp(firing_rate_contribution);

    disp('Variance Contributions (Factor 1, Factor 2):');
    disp(variance_contribution);

    % Compute dominance differences for each dataset
    firing_rate_dominance = abs(firing_rate_loadings(:, 1)) - abs(firing_rate_loadings(:, 2));
    variance_dominance = abs(variance_loadings(:, 1)) - abs(variance_loadings(:, 2));

    figure('Name','Dominance Differences (Factor 1 - Factor 2)','Color','w');
    subplot(1, 2, 1);
    plot(firing_rate_dominance, 'LineWidth',1.3);
    title('Firing Rate: Dominance Diff (F1 - F2)');
    xlabel('Time Slices');
    ylabel('Difference');

    subplot(1, 2, 2);
    plot(variance_dominance, 'LineWidth',1.3);
    title('Variance: Dominance Diff (F1 - F2)');
    xlabel('Time Slices');
    ylabel('Difference');

    % --------------------------------
    %  5) T-TESTS BETWEEN FACTOR LOADS
    % --------------------------------

    % Perform t-tests for Factor 1
    [h1, p1] = ttest(firing_rate_loadings(:, 1), variance_loadings(:, 1)); 
    % Perform t-tests for Factor 2
    [h2, p2] = ttest(firing_rate_loadings(:, 2), variance_loadings(:, 2)); 

    % Display t-test results
    disp('T-test results (Factor 1):');
    disp(['H = ', num2str(h1), ' (1 => significant difference)']);
    disp(['p-value = ', num2str(p1)]);

    disp('T-test results (Factor 2):');
    disp(['H = ', num2str(h2), ' (1 => significant difference)']);
    disp(['p-value = ', num2str(p2)]);

    % ------------------------
    %  6) BOXPLOTS OF LOADINGS
    % ------------------------

    % Boxplot for Factor 1
    figure('Name','Boxplots: Factor 1','Color','w');
    boxplot([firing_rate_loadings(:, 1), variance_loadings(:, 1)], ...
        'Labels', {'Firing Rate (F1)', 'Variance (F1)'});
    title('Comparison of Factor 1 Loadings');

    % Boxplot for Factor 2
    figure('Name','Boxplots: Factor 2','Color','w');
    boxplot([firing_rate_loadings(:, 2), variance_loadings(:, 2)], ...
        'Labels', {'Firing Rate (F2)', 'Variance (F2)'});
    title('Comparison of Factor 2 Loadings');

    % -----------------------------------
    %  7) CORRELATION BETWEEN FACTOR LOADS
    % -----------------------------------

    % Correlation between FR and variance loadings (Factor 1)
    correlation_factor1 = corr(firing_rate_loadings(:, 1), variance_loadings(:, 1));
    % Correlation between FR and variance loadings (Factor 2)
    correlation_factor2 = corr(firing_rate_loadings(:, 2), variance_loadings(:, 2));

    % Display correlation results
    disp('Correlation between Firing Rate and Variance Loadings:');
    disp(['Factor 1: ', num2str(correlation_factor1)]);
    disp(['Factor 2: ', num2str(correlation_factor2)]);

    % Scatter plot for Factor 1
    figure('Name','Correlation: Factor 1','Color','w');
    scatter(firing_rate_loadings(:, 1), variance_loadings(:, 1), 'filled');
    xlabel('Firing Rate Loadings (F1)');
    ylabel('Variance Loadings (F1)');
    title('Correlation: Factor 1');

    % Scatter plot for Factor 2
    figure('Name','Correlation: Factor 2','Color','w');
    scatter(firing_rate_loadings(:, 2), variance_loadings(:, 2), 'filled');
    xlabel('Firing Rate Loadings (F2)');
    ylabel('Variance Loadings (F2)');
    title('Correlation: Factor 2');

    % -------------------------------------------------------
    %  8) REGRESSION & RESIDUAL ANALYSIS (FIRING RATE -> VAR)
    % -------------------------------------------------------

    % Factor 1: Regress variance on firing rate
    [b,~,residuals] = regress(variance_loadings(:, 1), ...
                              [firing_rate_loadings(:, 1), ...
                               ones(size(firing_rate_loadings, 1), 1)]);
    disp('Regression Coefficients (Var Factor 1 ~ FR Factor 1):');
    disp(b);

    % Plot residuals
    figure('Name','Residuals: Factor 1','Color','w');
    plot(residuals, 'r','LineWidth',1.3);
    title('Residuals of Variance (Factor 1) After Removing FR Contribution');
    xlabel('Time Slices');
    ylabel('Residuals');

    % Factor 2: Regress variance on firing rate
    [b2,~,residuals2] = regress(variance_loadings(:, 2), ...
                                [firing_rate_loadings(:, 2), ...
                                 ones(size(firing_rate_loadings, 1), 1)]);
    figure('Name','Residuals: Factor 2','Color','w');
    plot(residuals2, 'b','LineWidth',1.3);
    title('Residuals of Variance (Factor 2) After Removing FR Contribution');
    xlabel('Time Slices');
    ylabel('Residuals');

    % --------------------------------
    %  9) PARTIAL CORRELATION ANALYSIS
    % --------------------------------

    % Partial correlation between Var(F1) and FR(F1), controlling for FR(F2)
    [partial_corr_factor1, pval_factor1] = partialcorr( ...
        variance_loadings(:, 1), ...
        firing_rate_loadings(:, 1), ...
        firing_rate_loadings(:, 2));

    % Partial correlation between Var(F2) and FR(F2), controlling for FR(F1)
    [partial_corr_factor2, pval_factor2] = partialcorr( ...
        variance_loadings(:, 2), ...
        firing_rate_loadings(:, 2), ...
        firing_rate_loadings(:, 1));

    disp('Partial Correlation (Var ~ FR) Controlling for Other Factor:');
    disp(['Factor 1: ', num2str(partial_corr_factor1), ...
          ', pval = ', num2str(pval_factor1)]);
    disp(['Factor 2: ', num2str(partial_corr_factor2), ...
          ', pval = ', num2str(pval_factor2)]);

    % ------------------------------------------------
    % 10) REPEAT FACTOR ANALYSIS WITH MORE FACTORS (3)
    % ------------------------------------------------

    num_factors = 3;
    [loadings_extra, psi_extra, stats_extra] = factoran(combined_data, num_factors);

    % Separate new loadings
    firing_rate_loadings_extra = loadings_extra(1:n_time_slices, :);
    variance_loadings_extra = loadings_extra(n_time_slices+1:end, :);

    figure('Name','FA Loadings (3 Factors)','Color','w');
    subplot(1, 2, 1);
    bar(firing_rate_loadings_extra);
    xlabel('Time Slices');
    ylabel('Loadings');
    title('Firing Rate Loadings (3 Factors)');

    subplot(1, 2, 2);
    bar(variance_loadings_extra);
    xlabel('Time Slices');
    ylabel('Loadings');
    title('Variance Loadings (3 Factors)');

    % -------------------------------------------------
    % 11) EXAMPLE NONLINEAR REGRESSION (ONE FACTOR->VAR)
    % -------------------------------------------------

    % Fit a nonlinear regression model, e.g. Var(F1) ~ a * FR(F1)^b + c
    nonlinear_model = fit(firing_rate_loadings(:, 1), ...
                          variance_loadings(:, 1), ...
                          'a*x^b + c');
    disp('Nonlinear Regression Model (Factor 1):');
    disp(nonlinear_model);

    figure('Name','Nonlinear Fit: Factor 1','Color','w');
    scatter(firing_rate_loadings(:, 1), variance_loadings(:, 1), 'filled');
    hold on;
    plot(nonlinear_model, 'r-');
    xlabel('Firing Rate Loadings (Factor 1)');
    ylabel('Variance Loadings (Factor 1)');
    title('Nonlinear Relationship (Factor 1)');
    hold off;

    % Factor 2
    nonlinear_model2 = fit(firing_rate_loadings(:, 2), ...
                           variance_loadings(:, 2), ...
                           'a*x^b + c');
    disp('Nonlinear Regression Model (Factor 2):');
    disp(nonlinear_model2);

    figure('Name','Nonlinear Fit: Factor 2','Color','w');
    scatter(firing_rate_loadings(:, 2), variance_loadings(:, 2), 'filled');
    hold on;
    plot(nonlinear_model2, 'r-');
    xlabel('Firing Rate Loadings (Factor 2)');
    ylabel('Variance Loadings (Factor 2)');
    title('Nonlinear Relationship (Factor 2)');
    hold off;

    % ---------------------------
    % 12) PCA ON THE COMBINED DATA
    % ---------------------------

    % PCA for firing-rate columns only
    [coeff_fr, score_fr, explained_fr] = pca(combined_data(:, 1:n_time_slices));
    % PCA for variance columns only
    [coeff_v, score_v, explained_v] = pca(combined_data(:, n_time_slices+1:end));

    % Visualize the variance explained by each component
    figure('Name','PCA Explained Variance','Color','w');
    subplot(2,1,1)
    bar(explained_fr);
    xlabel('Principal Components');
    ylabel('Firing Rate Explained (%)');
    title('PCA: Firing Rate Explained by Components');

    subplot(2,1,2)
    bar(explained_v);
    xlabel('Principal Components');
    ylabel('Variance Explained (%)');
    title('PCA: Variance Explained by Components');

    % Scatter plot of the first three principal component scores
    figure('Name','PCA 3D Scores','Color','w');
    scatter3(score_fr(:, 1), score_fr(:, 2), score_fr(:, 3), 'filled');
    hold on;
    scatter3(score_v(:, 1), score_v(:, 2), score_v(:, 3), 'filled');
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    title('PCA: First Three Principal Components (FR & Var)');
    legend({'Firing Rate','Variance'});
    hold off;

end
