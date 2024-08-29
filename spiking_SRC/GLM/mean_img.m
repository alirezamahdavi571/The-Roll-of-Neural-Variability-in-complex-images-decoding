function mean = mean_img(img) 
    % Flatten the image to a single vector
    img_flatten = img(:);
    sum_t = sum(img_flatten);
    mean = sum_t/length(img_flatten);
end