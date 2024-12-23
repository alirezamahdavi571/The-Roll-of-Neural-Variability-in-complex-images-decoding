function [accuracy_w, accuracy_w_mat] = lda_classifier_var_mean(var_vec, number_of_time_slices, region)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SVM Classifier Function on Variance and Mean
    % Inputs:
    %   - var_vec: Input feature matrix of dimensions (samples, features, time slices)
    %   - number_of_time_slices: Total number of time slices
    %   - region: String indicating the brain region or data region
    % Outputs:
    %   - accuracy_w: Accuracy for each time slice
    %   - recall_artifact: Recall for artifact features across time slices
    %   - recall_body: Recall for body features across time slices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    iteration = 10; % Number of iterations

    % Extract features for each class
    fanofactor_artifact_features = var_vec(:, 1:36, :);
    fanofactor_body_features = var_vec(:, 37:72, :);
    fanofactor_natural_features = var_vec(:, 73:108, :);
    fanofactor_face_features = var_vec(:, 109:144, :);

    % Initialize outputs
    out = struct();
    accuracy_w = zeros(number_of_time_slices, 1);

    % SVM classification for each time slice
    for ts = 1:number_of_time_slices
        % Prepare features
        artifact_features = transpose(fanofactor_artifact_features(3:end, :, ts));
        body_features = transpose(fanofactor_body_features(3:end, :, ts));
        natural_features = transpose(fanofactor_natural_features(3:end, :, ts));
        face_features = transpose(fanofactor_face_features(3:end, :, ts));

        % Create labels
        artifact_labels = -1 * ones(size(artifact_features, 1), 1);
        body_labels = 1 * ones(size(body_features, 1), 1);
        natural_labels = -1 * ones(size(natural_features, 1), 1);
        face_labels = 1 * ones(size(face_features, 1), 1);

        % Stack features and labels
        features = [artifact_features; body_features; natural_features; face_features];
        labels = [artifact_labels; body_labels; natural_labels; face_labels];

        % Perform SVM decoding
        disp(['Accuracy for time slice ' num2str(ts)]);
        out(ts).out = gen_fx_get_lda(labels, features, 0.5, iteration);
        accuracy_w(ts) = mean(out(ts).out.pt, 1) * 100;
        accuracy_w_mat(ts,:) = out(ts).out.pt * 100;
    end

    % Time vector for plotting
    time = linspace(-200, 700, number_of_time_slices);

    % Plot accuracy
    figure();
    plot(time, accuracy_w, '-o', 'LineWidth', 5);
    hold on;
    plot(time, movmean(accuracy_w, 10), 'LineWidth', 5);
    grid on;
    legend({'Accuracy', 'Smoothed Accuracy'}, 'FontSize', 16, 'Location', 'best');
    legend('boxoff');
    [t, s] = title(['Accuracy on Mean in ', region], 'Color', 'black');
    t.FontSize = 16;
    s.FontAngle = 'italic';
    xlabel('Time (ms)', 'FontSize', 16, 'Color', 'b');
    ylabel('Accuracy Percentage', 'FontSize', 16, 'Color', 'r');
    save(['Accuracy_animate_inanimate_var_', region, '.mat'], 'accuracy_w', '-v7.3');

    % Compute recall
    recall_artifact = zeros(1, number_of_time_slices);
    recall_body = zeros(1, number_of_time_slices);
    for ts = 1:number_of_time_slices
        confusion_mat = out(ts).out.C;
        recall_artifact_temp = zeros(1, iteration);
        recall_body_temp = zeros(1, iteration);
        for i = 1:iteration
            temp_confusion_matrix = confusion_mat(:, :, i);
            recall_artifact_temp(i) = temp_confusion_matrix(1, 1) / (temp_confusion_matrix(1, 1) + temp_confusion_matrix(1, 2));
            recall_body_temp(i) = temp_confusion_matrix(2, 2) / (temp_confusion_matrix(2, 1) + temp_confusion_matrix(2, 2));
        end
        recall_artifact(ts) = mean(recall_artifact_temp, 2) * 100;
        recall_body(ts) = mean(recall_body_temp, 2) * 100;
    end

    % Plot recall
    figure();
    plot(time, movmean(recall_artifact, 10), 'LineWidth', 5);
    hold on;
    plot(time, movmean(recall_body, 10), 'LineWidth', 5);
    grid on;
    legend({'Recall Inanimate', 'Recall Animate'}, 'FontSize', 16, 'Location', 'best');
    legend('boxoff');
    [t, s] = title(['Recall on Variance in ', region], 'Color', 'black');
    t.FontSize = 16;
    s.FontAngle = 'italic';
    xlabel('Time (ms)', 'FontSize', 16, 'Color', 'b');
    ylabel('Recall Percentage', 'FontSize', 16, 'Color', 'r');
    save(['Recall_animate_var_', region, '.mat'], 'recall_artifact', '-v7.3');
    save(['Recall_inanimate_var_', region, '.mat'], 'recall_body', '-v7.3');
end
