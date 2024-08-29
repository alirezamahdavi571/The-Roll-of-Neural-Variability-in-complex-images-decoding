function SpikeTrain_it_all = seperateStimulus(SpikeTrain_it_all,catlabels,catagory,number_of_neurons)
    for i = 3:number_of_neurons
       counter = 0;
       cm =  SpikeTrain_it_all(i).cm;
       
       
       switch(catagory)
           case 'face'
                for j = 1:length(cm)
                    for k = 1:length(catlabels)
                        if(catlabels(k) == cm(j))
                            counter = counter +1;
                            SpikeTrain_it_all(i).cm_face(counter) = cm(j);
                            SpikeTrain_it_all(i).cm_face_idx(counter) = j ;
                        end
                    end
                end
           case 'body'
                for j = 1:length(cm)
                    for k = 1:length(catlabels)
                        if(catlabels(k) == cm(j))
                            counter = counter +1;
                            SpikeTrain_it_all(i).cm_body(counter) = cm(j);
                            SpikeTrain_it_all(i).cm_body_idx(counter) = j ;
                        end
                    end
                end
           case 'natural'
                for j = 1:length(cm)
                    for k = 1:length(catlabels)
                        if(catlabels(k) == cm(j))
                            counter = counter +1;
                            SpikeTrain_it_all(i).cm_natural(counter) = cm(j);
                            SpikeTrain_it_all(i).cm_natural_idx(counter) = j ;
                        end
                    end
                end
   
           case 'artifact'
                for j = 1:length(cm)
                    for k = 1:length(catlabels)
                        if(catlabels(k) == cm(j))
                            counter = counter +1;
                            SpikeTrain_it_all(i).cm_artfact(counter) = cm(j);
                            SpikeTrain_it_all(i).cm_artifact_idx(counter) = j ;
                        end
                    end
                end
           case 'nonface'
                 for j = 1:length(cm)
                     for k = 1:length(catlabels)
                        if(catlabels(k) == cm(j))
                            counter = counter +1;
                            SpikeTrain_it_all(i).cm_nonface(counter) = cm(j);
                            SpikeTrain_it_all(i).cm_nonface_idx(counter) = j ;
                        end
                     end
                end
           otherwise
               disp('category not found!!!');
       end
    end
end