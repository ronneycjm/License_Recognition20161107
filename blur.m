I = imread('test_blur.jpg');
PSF = fspecial('motion',20, 30);
im_blur = imfilter(I, PSF, 'circular','conv');
imwrite(im_blur,'test_blur.jpg');

