function [accuracy_combined,accuracy_combined_matrix] = joint_lda_classifier(var_vec, artifact_Mean_PSTH_data, body_Mean_PSTH_data, natural_Mean_PSTH_data, Face_Mean_PSTH_data, number_of_time_slices, iteration)
    % Function to compute feature-based SVM decoding for each time slice.
    %
    % Inputs:
    %   var_vec - Variance (fanofactor) data matrix
    %   artifact_Mean_PSTH_data, body_Mean_PSTH_data, natural_Mean_PSTH_data, Face_Mean_PSTH_data - Firing rate data matrices
    %   number_of_time_slices - Total number of time slices
    %   iteration - Number of iterations for LDA
    %
    % Outputs:
    %   accuracy_combined - Combined accuracy for each time slice

    % Initialize outputs
    accuracy_combined = zeros(number_of_time_slices, 1);
    out = struct();

    % Extract features for each class
    fanofactor_artifact_features = var_vec(:, 1:36, :);
    fanofactor_body_features = var_vec(:, 37:72, :);
    fanofactor_natural_features = var_vec(:, 73:108, :);
    fanofactor_face_features = var_vec(:, 109:144, :);

    % SVM Decoding for each time slice
    for ts = 1:number_of_time_slices
        % Extract features
        [artifact_features, body_features, natural_features, face_features] = extractFeatures(ts, ...
            artifact_Mean_PSTH_data, body_Mean_PSTH_data, natural_Mean_PSTH_data, Face_Mean_PSTH_data, ...
            fanofactor_artifact_features, fanofactor_body_features, fanofactor_natural_features, fanofactor_face_features);

        % Create labels
        [features, labels] = stackFeaturesAndLabels(artifact_features, body_features, natural_features, face_features);

        % Train SVM and compute accuracy
        disp(['Accuracy for time slice ' num2str(ts)]);
        out(ts).out = gen_fx_get_lda(labels, features, 0.5, iteration);
    end

    % Extract accuracies
    for ts = 1:number_of_time_slices
        accuracy_combined(ts, 1) = mean(out(ts).out.pt, 1) * 100;
        accuracy_combined_matrix(ts,:) = out(ts).out.pt * 100;
    end

    % Display final accuracy
    disp('Combined Features Accuracy per Time Slice:');
    disp(accuracy_combined);

    % Save the results
    save('svm_accuracy_per_time_slice.mat', 'accuracy_combined', '-v7.3');
end

function [artifact_features, body_features, natural_features, face_features] = extractFeatures(ts, artifact_Mean_PSTH_data, body_Mean_PSTH_data, natural_Mean_PSTH_data, Face_Mean_PSTH_data, fanofactor_artifact_features, fanofactor_body_features, fanofactor_natural_features, fanofactor_face_features)
    % Helper function to extract features for a specific time slice
    firing_rate_artifact = transpose(artifact_Mean_PSTH_data(3:end, :, ts));
    firing_rate_body = transpose(body_Mean_PSTH_data(3:end, :, ts));
    firing_rate_natural = transpose(natural_Mean_PSTH_data(3:end, :, ts));
    firing_rate_face = transpose(Face_Mean_PSTH_data(3:end, :, ts));
    
    var_artifact = transpose(fanofactor_artifact_features(3:end, :, ts));
    var_body = transpose(fanofactor_body_features(3:end, :, ts));
    var_natural = transpose(fanofactor_natural_features(3:end, :, ts));
    var_face = transpose(fanofactor_face_features(3:end, :, ts));
    
    artifact_features = [firing_rate_artifact, var_artifact];
    body_features = [firing_rate_body, var_body];
    natural_features = [firing_rate_natural, var_natural];
    face_features = [firing_rate_face, var_face];
end

function [features, labels] = stackFeaturesAndLabels(artifact_features, body_features, natural_features, face_features)
    % Helper function to stack features and create labels
    artifact_labels = 1 * ones(size(artifact_features, 1), 1);
    body_labels = -1 * ones(size(body_features, 1), 1);
    natural_labels = 1 * ones(size(natural_features, 1), 1);
    face_labels = -1 * ones(size(face_features, 1), 1);
    
    features = [artifact_features; body_features; natural_features; face_features];
    labels = [artifact_labels; body_labels; natural_labels; face_labels];
end
