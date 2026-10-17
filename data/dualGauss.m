%clear; close all; clc;

imageSize = 1024;  
k=imageSize/1024;
center = [imageSize/2, imageSize/2];
[X, Y] = meshgrid(1:imageSize, 1:imageSize);

mu_g1 = center - [70, -50]*3*k;
mu_g2 = center - [-80, 60]*3*k;
sigma_g1 = [20000, 8000; 8000, 25000]*2*k^2;
sigma_g2 = [40000, -12000; -12000, 30000]*2*k^2;
 
x_g1 = X - mu_g1(1);
y_g1 = Y - mu_g1(2);
sigma_inv1 = inv(sigma_g1);
Z_g1 = exp(-0.5 * (sigma_inv1(1,1)*x_g1.^2 + 2*sigma_inv1(1,2)*x_g1.*y_g1 + sigma_inv1(2,2)*y_g1.^2));
 
x_g2 = X - mu_g2(1);
y_g2 = Y - mu_g2(2);
sigma_inv2 = inv(sigma_g2);
Z_g2 = exp(-0.5 * (sigma_inv2(1,1)*x_g2.^2 + 2*sigma_inv2(1,2)*x_g2.*y_g2 + sigma_inv2(2,2)*y_g2.^2));
 
Z = (Z_g1 + Z_g2)/2;
Z = Z / max(Z(:));  
cauchyImage = uint8(Z * 255);
Z=Z/(sum(sum(Z)));
filename = 'dualgauss.mat';
save(filename, 'Z');
figure('Position', [100, 100, 800, 800]);
imshow(cauchyImage);
%title('Test Measure: 2D Cauchy Distribution (1024x1024)', 'FontSize', 14);
colormap('parula');  
colorbar;
% imwrite(cauchyImage, 'cauchy_1024.png');