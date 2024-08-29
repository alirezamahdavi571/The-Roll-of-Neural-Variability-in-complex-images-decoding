function contrast = contrast(img) 

    % Flatten the image to a single vector
    img_flatten = img(:);

    % Sort the pixel values
    sorted_pixels = sort(img_flatten);

    % Find the quartile indices
    n = length(sorted_pixels);
    quartile_index = floor(n / 4);

    % Calculate the first and third quartiles
    qrt1 = sorted_pixels(1:quartile_index+1);
    qrt3 = sorted_pixels(3*quartile_index:end);

    % Calculate the custom contrast measure
    contrast = (mean(qrt3) - mean(qrt1)) / (mean(qrt3) + mean(qrt1));
end