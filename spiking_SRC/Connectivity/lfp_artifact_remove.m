function artif=lfp_artifact_remove(sig,preview)

sigma=1.5;
artif=[];
if isempty(preview)
preview=0;
end


    
    
   
        val_it_max=max(sig,[],2);
        val_it_min=min(sig,[],2);
        val_it=val_it_max-val_it_min;
        th_it=median(val_it)+sigma*std(val_it);
        val_it=max(sig,[],2)-min(sig,[],2);
        artif=zeros(size(sig,1),1);
        artif(abs(val_it)>th_it)=1;
        artif(abs(val_it)==0)=1; 
        
        artif(abs(val_it_max)>th_it)=1;
        artif(abs(val_it_min)>th_it)=1;

   
    
       
    
    
    if preview
        figure
        subplot(211)
        plot((sig(artif==1,:))')
         
         subplot(212)
        plot((sig(artif==0,:))')
      
    end
    
    
end


