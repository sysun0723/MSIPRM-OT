%clear; close all; clc;

imageSize = 1024;  
k=imageSize/1024;
center = [imageSize/2, imageSize/2];
[X, Y] = meshgrid(1:imageSize, 1:imageSize);

mu_g = center+[100,-100]*k;
sigma_g = [30000, 10000; 20000, 30000]*4*k^2; 
 
x_g = X - mu_g(1);
y_g = Y - mu_g(2);
sigma_inv = inv(sigma_g);
Z = exp(-0.5 * (sigma_inv(1,1)*x_g.^2 + 2*sigma_inv(1,2)*x_g.*y_g + sigma_inv(2,2)*y_g.^2));
Z = Z / max(Z(:));  
gaussImage = uint8(Z * 255);

Z=Z/(sum(sum(Z)));

filename = 'Singlegauss.mat';
save(filename, 'Z');
figure('Position', [100, 100, 800, 800]);
imshow(gaussImage);
%title('Test Measure: 2D Cauchy Distribution (1024x1024)', 'FontSize', 14);
colormap('parula');  
colorbar;
% imwrite(cauchyImage, 'cauchy_1024.png');

