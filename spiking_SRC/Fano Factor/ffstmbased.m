function [mean_vec , var_vec] = ffstmbased(SpikeTrain_it_all , mean_vec,var_vec,All_stimuli,number_of_neurons,number_of_time_slices,sliding_step,window_length)
    for i = 3:number_of_neurons
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
        cm =  SpikeTrain_it_all(i).cm;
        data = SpikeTrain_it_all(i).data;
        for j = 1:length(All_stimuli)
            test = All_stimuli(j) == cm;
            data_undertest_per_neroun = data(test,:);
         
            for k = 1:size(data_undertest_per_neroun,1)
                for u = 1:number_of_time_slices
                    temp = data_undertest_per_neroun(k,1+sliding_step*(u-1):window_length+sliding_step*(u-1));
                    spike_count_stimuli(k,u) = sum(temp);
                end
            end
        
            for u = 1:number_of_time_slices
                mean_vec(i,j,u) = mean(spike_count_stimuli(:,u));
                var_vec(i,j,u)  = var(spike_count_stimuli(:,u));
            end
        end
   
   
    end
end