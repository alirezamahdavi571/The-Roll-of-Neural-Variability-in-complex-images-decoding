function tt_decoding_var = time_time_decoding_VAR(var_vec, number_of_time_slices, region)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time-Time Decoding Function
    % Inputs:
    %   - var_vec: Input feature matrix of dimensions (samples, features, time slices)
    %   - number_of_time_slices: Total number of time slices
    %   - region: String indicating the brain region or data region
    % Output:
    %   - tt_decoding_var: Decoding matrix for all time slices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Extract features for each class
    fanofactor_all_features = var_vec;
    fanofactor_artifact_features = fanofactor_all_features(:, 1:36, :);
    fanofactor_body_features = fanofactor_all_features(:, 37:72, :);
    fanofactor_natural_features = fanofactor_all_features(:, 73:108, :);
    fanofactor_face_features = fanofactor_all_features(:, 109:144, :);
    
    % Initialize output structure
    iteration = 1;
    out = struct();
    
    % Decoding within each time slice
    for ts = 1:number_of_time_slices
        % Prepare features
        artifact_features = transpose(fanofactor_artifact_features(3:end, :, ts));
        body_features = transpose(fanofactor_body_features(3:end, :, ts));
        natural_features = transpose(fanofactor_natural_features(3:end, :, ts));
        face_features = transpose(fanofactor_face_features(3:end, :, ts));
        
        % Create labels
        artifact_labels = -1*ones(size(artifact_features, 1), 1);
        body_labels = -1 * ones(size(body_features, 1), 1);
        natural_labels = -1 * ones(size(natural_features, 1), 1);
        face_labels = 1 * ones(size(face_features, 1), 1);
        
        % Combine features and labels
        features = [artifact_features; body_features; natural_features; face_features];
        labels = [artifact_labels; body_labels; natural_labels; face_labels];
        
        % Perform decoding for current time slice
        disp(['Accuracy for time slice ' num2str(ts)]);
        out(ts).out = gen_fx_get_svm(labels, features, 0.5, iteration);
    end
    
    % Generalization across time slices
    out_generalization = struct();
    for j = 1:number_of_time_slices
        for ts = 1:number_of_time_slices
            % Prepare features
            artifact_features = transpose(fanofactor_artifact_features(3:end, :, ts));
            body_features = transpose(fanofactor_body_features(3:end, :, ts));
            natural_features = transpose(fanofactor_natural_features(3:end, :, ts));
            face_features = transpose(fanofactor_face_features(3:end, :, ts));
            
            % Create labels
            artifact_labels = -1*ones(size(artifact_features, 1), 1);
            body_labels = -1 * ones(size(body_features, 1), 1);
            natural_labels = -1 * ones(size(natural_features, 1), 1);
            face_labels = 1 * ones(size(face_features, 1), 1);
            
            % Combine features and labels
            features = [artifact_features; body_features; natural_features; face_features];
            labels = [artifact_labels; body_labels; natural_labels; face_labels];
            
            % Perform generalization decoding
            disp(['Accuracy for time slice ' num2str(j) ' on ts ' num2str(ts)]);
            out_generalization(j, ts).out = gen_fx_get_svm_EEG(labels, features, 0.5, iteration, 'generalization', out(j).out.model);
        end
    end
    
    % Collect decoding results into a matrix
    tt_decoding_var = zeros(number_of_time_slices, number_of_time_slices);
    for j = 1:number_of_time_slices
        for ts = 1:number_of_time_slices
            tt_decoding_var(j, ts) = out_generalization(j, ts).out.pt;
        end
    end
    
    % Save results
    save(['time_time_decoding_var_face_', region, '.mat'], 'tt_decoding_var', '-v7.3');
end
