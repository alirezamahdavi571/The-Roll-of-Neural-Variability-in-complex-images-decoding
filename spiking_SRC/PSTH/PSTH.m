function SpikeTrain_it_all = PSTH(SpikeTrain_it_all,number_of_neurons,category,sliding_step,window_length,number_of_time_slices)
    dt = 0.001;
    for i = 3:number_of_neurons
    
        switch(category)
            case 'face'
                SpikeTrain_it_face = SpikeTrain_it_all(i).data;
                cm_face =  unique(SpikeTrain_it_all(i).cm_face);
                cm = SpikeTrain_it_all(i).cm;
 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %First step
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for j = 1:length(cm_face)
                    trials_face = find(cm == cm_face(j));
                    raster_face = zeros(length(trials_face) , 900);
        
                    for m = 1:length(trials_face)
                        raster_face(m,:) = SpikeTrain_it_face(trials_face(m),:);
                    end
                
                    SpikeTrain_it_all(i).rasters_face(j).data = raster_face;
  
    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Second step
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    raster_Under_train_face = SpikeTrain_it_all(i).rasters_face(j).data;
                    for m = 1:size(raster_Under_train_face,2)
                        counter_face = 0;
                        for k = 1:size(raster_Under_train_face,1)
                            if(raster_Under_train_face(k,m) == 1)
                                counter_face = counter_face + 1 ;
                            end
                        end
                        PSTH_face(j,m) = counter_face /((dt)*size(raster_Under_train_face,1));
                    end
                    PSTH_face(j,:) = movmean(PSTH_face(j,:) , 100);
                end

                SpikeTrain_it_all(i).psth_face = PSTH_face;
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %third step
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for j = 1:length(cm_face)
                    for u = 1: number_of_time_slices
                        temp_face = SpikeTrain_it_all(i).psth_face(j,1+sliding_step*(u-1):window_length+sliding_step*(u-1));
                        Face_Mean_PSTH_data(j ,u) = mean(temp_face);
                        Face_Var_PSTH_data(j,u)= var(temp_face);
                    end
                end
    
                SpikeTrain_it_all(i).Face_Mean_PSTH_data = Face_Mean_PSTH_data;
                SpikeTrain_it_all(i).Face_Var_PSTH_data = Face_Var_PSTH_data;
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % windowing and moving mean on PSTH per Stimuli per window any window has
                % one SVM classifier 
                % shape is = (number_of_neurons , number_of_Stimuli , number_of_time_slices)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            case 'natural'
                SpikeTrain_it_natural = SpikeTrain_it_all(i).data;
                cm_natural =  unique(SpikeTrain_it_all(i).cm_natural);
                cm = SpikeTrain_it_all(i).cm;
 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %First step
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for j = 1:length(cm_natural)
                    trials_natural = find(cm == cm_natural(j));
                    raster_natural = zeros(length(trials_natural) , 900);
        
                    for m = 1:length(trials_natural)
                        raster_natural(m,:) = SpikeTrain_it_natural(trials_natural(m),:);
                    end
                
                    SpikeTrain_it_all(i).rasters_natural(j).data = raster_natural;
  
    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Second step
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    raster_Under_train_natural = SpikeTrain_it_all(i).rasters_natural(j).data;
                    for m = 1:size(raster_Under_train_natural,2)
                        counter_natural = 0;
                        for k = 1:size(raster_Under_train_natural,1)
                            if(raster_Under_train_natural(k,m) == 1)
                                counter_natural = counter_natural + 1 ;
                            end
                        end
                        PSTH_natural(j,m) = counter_natural /((dt)*size(raster_Under_train_natural,1));
                    end
                    PSTH_natural(j,:) = movmean(PSTH_natural(j,:) , 100);
                end

                SpikeTrain_it_all(i).psth_natural = PSTH_natural;
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %third step
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for j = 1:length(cm_natural)
                    for u = 1: number_of_time_slices
                        temp_natural = SpikeTrain_it_all(i).psth_natural(j,1+sliding_step*(u-1):window_length+sliding_step*(u-1));
                        natural_Mean_PSTH_data(j ,u) = mean(temp_natural);
                        natural_Var_PSTH_data(j,u)= var(temp_natural);
                    end
                end
    
                SpikeTrain_it_all(i).natural_Mean_PSTH_data = natural_Mean_PSTH_data;
                SpikeTrain_it_all(i).natural_Var_PSTH_data = natural_Var_PSTH_data;
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % windowing and moving mean on PSTH per Stimuli per window any window has
                % one SVM classifier 
                % shape is = (number_of_neurons , number_of_Stimuli , number_of_time_slices)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'body'
                SpikeTrain_it_body = SpikeTrain_it_all(i).data;
                cm_body =  unique(SpikeTrain_it_all(i).cm_body);
                cm = SpikeTrain_it_all(i).cm;
 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %First step
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for j = 1:length(cm_body)
                    trials_body = find(cm == cm_body(j));
                    raster_body = zeros(length(trials_body) , 900);
        
                    for m = 1:length(trials_body)
                        raster_body(m,:) = SpikeTrain_it_body(trials_body(m),:);
                    end
                
                    SpikeTrain_it_all(i).rasters_body(j).data = raster_body;
  
    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Second step
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    raster_Under_train_body = SpikeTrain_it_all(i).rasters_body(j).data;
                    for m = 1:size(raster_Under_train_body,2)
                        counter_body = 0;
                        for k = 1:size(raster_Under_train_body,1)
                            if(raster_Under_train_body(k,m) == 1)
                                counter_body = counter_body + 1 ;
                            end
                        end
                        PSTH_body(j,m) = counter_body /((dt)*size(raster_Under_train_body,1));
                    end
                    PSTH_body(j,:) = movmean(PSTH_body(j,:) , 100);
                end

                SpikeTrain_it_all(i).psth_body = PSTH_body;
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %third step
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for j = 1:length(cm_body)
                    for u = 1: number_of_time_slices
                        temp_body = SpikeTrain_it_all(i).psth_body(j,1+sliding_step*(u-1):window_length+sliding_step*(u-1));
                        body_Mean_PSTH_data(j ,u) = mean(temp_body);
                        body_Var_PSTH_data(j,u)= var(temp_body);
                    end
                end
    
                SpikeTrain_it_all(i).body_Mean_PSTH_data = body_Mean_PSTH_data;
                SpikeTrain_it_all(i).body_Var_PSTH_data = body_Var_PSTH_data;
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % windowing and moving mean on PSTH per Stimuli per window any window has
                % one SVM classifier 
                % shape is = (number_of_neurons , number_of_Stimuli , number_of_time_slices)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'artifact'
                SpikeTrain_it_artifact = SpikeTrain_it_all(i).data;
                cm_artifact =  unique(SpikeTrain_it_all(i).cm_artifact);
                cm = SpikeTrain_it_all(i).cm;
 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %First step
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for j = 1:length(cm_artifact)
                    trials_artifact = find(cm == cm_artifact(j));
                    raster_artifact = zeros(length(trials_artifact) , 900);
        
                    for m = 1:length(trials_artifact)
                        raster_artifact(m,:) = SpikeTrain_it_artifact(trials_artifact(m),:);
                    end
                
                    SpikeTrain_it_all(i).rasters_artifact(j).data = raster_artifact;
  
    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Second step
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    raster_Under_train_artifact = SpikeTrain_it_all(i).rasters_artifact(j).data;
                    for m = 1:size(raster_Under_train_artifact,2)
                        counter_artifact = 0;
                        for k = 1:size(raster_Under_train_artifact,1)
                            if(raster_Under_train_artifact(k,m) == 1)
                                counter_artifact = counter_artifact + 1 ;
                            end
                        end
                        PSTH_artifact(j,m) = counter_artifact /((dt)*size(raster_Under_train_artifact,1));
                    end
                    PSTH_artifact(j,:) = movmean(PSTH_artifact(j,:) , 100);
                end

                SpikeTrain_it_all(i).psth_artifact = PSTH_artifact;
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %third step
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for j = 1:length(cm_artifact)
                    for u = 1: number_of_time_slices
                        temp_artifact = SpikeTrain_it_all(i).psth_artifact(j,1+sliding_step*(u-1):window_length+sliding_step*(u-1));
                        artifact_Mean_PSTH_data(j ,u) = mean(temp_artifact);
                        artifact_Var_PSTH_data(j,u)= var(temp_artifact);
                    end
                end
    
                SpikeTrain_it_all(i).artifact_Mean_PSTH_data = artifact_Mean_PSTH_data;
                SpikeTrain_it_all(i).artifact_Var_PSTH_data = artifact_Var_PSTH_data;
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % windowing and moving mean on PSTH per Stimuli per window any window has
                % one SVM classifier 
                % shape is = (number_of_neurons , number_of_Stimuli , number_of_time_slices)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'nonface'
                SpikeTrain_it_nonface = SpikeTrain_it_all(i).data;
                cm_nonface =  unique(SpikeTrain_it_all(i).cm_nonface);
                cm = SpikeTrain_it_all(i).cm;
 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %First step
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for j = 1:length(cm_nonface)
                    trials_nonface = find(cm == cm_nonface(j));
                    raster_nonface = zeros(length(trials_nonface) , 900);
        
                    for m = 1:length(trials_nonface)
                        raster_nonface(m,:) = SpikeTrain_it_nonface(trials_nonface(m),:);
                    end
                
                    SpikeTrain_it_all(i).rasters_nonface(j).data = raster_nonface;
  
    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Second step
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    raster_Under_train_nonface = SpikeTrain_it_all(i).rasters_nonface(j).data;
                    for m = 1:size(raster_Under_train_nonface,2)
                        counter_nonface = 0;
                        for k = 1:size(raster_Under_train_nonface,1)
                            if(raster_Under_train_nonface(k,m) == 1)
                                counter_nonface = counter_nonface + 1 ;
                            end
                        end
                        PSTH_nonface(j,m) = counter_nonface /((dt)*size(raster_Under_train_nonface,1));
                    end
                    PSTH_nonface(j,:) = movmean(PSTH_nonface(j,:) , 100);
                end

                SpikeTrain_it_all(i).psth_nonface = PSTH_nonface;
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %third step
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for j = 1:length(cm_nonface)
                    for u = 1: number_of_time_slices
                        temp_nonface = SpikeTrain_it_all(i).psth_nonface(j,1+sliding_step*(u-1):window_length+sliding_step*(u-1));
                        nonface_Mean_PSTH_data(j ,u) = mean(temp_nonface);
                        nonface_Var_PSTH_data(j,u)= var(temp_nonface);
                    end
                end
    
                SpikeTrain_it_all(i).nonface_Mean_PSTH_data = nonface_Mean_PSTH_data;
                SpikeTrain_it_all(i).nonface_Var_PSTH_data = nonface_Var_PSTH_data;
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % windowing and moving mean on PSTH per Stimuli per window any window has
                % one SVM classifier 
                % shape is = (number_of_neurons , number_of_Stimuli , number_of_time_slices)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            otherwise
                disp('catagory not found!!!');
        end

    end

end
