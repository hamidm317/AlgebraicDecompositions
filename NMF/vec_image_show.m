function Done = vec_image_show(vec_image, m, n)

    img = uint8(255 * reshape(vec_image, m, n));
    
    imshow(img)
    
    Done = 1;

end

