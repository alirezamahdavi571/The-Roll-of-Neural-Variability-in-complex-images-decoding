function [spikeItAll,fanofactor] = fanoFactorAll(SpikeTrain_it_all,number_of_neurons,number_of_time_slices,window_length,sliding_step)

    mean_vec = zeros(number_of_neurons,number_of_time_slices);
    var_vec = zeros(number_of_neurons,number_of_time_slices);
    spike_count = zeros(number_of_neurons,900);
    temp = zeros(1,window_length);

    figure()
    spy(SpikeTrain_it_all(20).data);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIND MEAN & VARIANCE FOR ANY NEURONS (WE HAVE 859 NEURONS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 3:number_of_neurons
        SpikeTrain_it = SpikeTrain_it_all(i).data;
        SpikeTrain_it_all(i).all_Trial_size = size(SpikeTrain_it,1);
        spike_count_2 = zeros(size(SpikeTrain_it,1),number_of_time_slices);
        spy(SpikeTrain_it);
    
        counter_face = 0 ; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COUNTING NUMBER OF SPIKES 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        for j = 1:size(SpikeTrain_it,1)
            for u = 1:number_of_time_slices
                temp = SpikeTrain_it(j,1+sliding_step*(u-1):window_length+sliding_step*(u-1));
                spike_count_2(j,u) = sum(temp);
            end
        end
    
    
        SpikeTrain_it_all(i).countSpikes = spike_count_2;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE MEAN AND VAR 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        for u = 1:number_of_time_slices
            mean_vec(i,u) = mean(spike_count_2(:,u));
            var_vec(i,u)  = var(spike_count_2(:,u));
        end
    end

    for i = 1:number_of_time_slices
        [fanofactor(i),~,~] = regression(mean_vec(:,i), var_vec(:,i),'one');
        if( i == 1)
            figure()
            subplot(2,1,1)
            scatter(mean_vec(:,i) , var_vec(:,i) , '*');title('time slice 1');xlabel('Mean');ylabel('Variance')
            grid on

            subplot(2,1,2)
            histogram(mean_vec(:,i));title('mean Dist')
            grid on
        end
    
    
        if( i == 20)
            figure()
            subplot(2,1,1)
            scatter(mean_vec(:,i) , var_vec(:,i) , '*');title('time slice 20');xlabel('Mean');ylabel('Variance')
            grid on

            subplot(2,1,2)
            histogram(mean_vec(:,i));title('mean Dist')
            grid on
        end


    end
    
    spikeItAll = SpikeTrain_it_all;
end