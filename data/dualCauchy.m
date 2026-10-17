%clear; close all; clc;
 
imageSize = 1024;  
k=imageSize/1024;
center = [imageSize/2, imageSize/2];
[X, Y] = meshgrid(1:imageSize, 1:imageSize);

mu_c1 = center - [400, 400]*0.75*k;
mu_c2 = center + [400, 400]*0.75*k;
gamma_c1 = [500, 400]*k;
gamma_c2 = [500, 600]*k;
rho_c1 = 0.4;
rho_c2 = -0.3;

x_c1 = (X - mu_c1(1)) / gamma_c1(1);
y_c1 = (Y - mu_c1(2)) / gamma_c1(2);
Z_c1 = 1 ./ (2*pi*sqrt(1-rho_c1^2)) .* ...
       1 ./ (1 + (x_c1.^2 - 2*rho_c1*x_c1.*y_c1 + y_c1.^2)/(1-rho_c1^2)).^(3/2);

x_c2 = (X - mu_c2(1)) / gamma_c2(1);
y_c2 = (Y - mu_c2(2)) / gamma_c2(2);
Z_c2 = 1 ./ (2*pi*sqrt(1-rho_c2^2)) .* ...
       1 ./ (1 + (x_c2.^2 - 2*rho_c2*x_c2.*y_c2 + y_c2.^2)/(1-rho_c2^2)).^(3/2);
 
Z = (Z_c1 + Z_c2)/2;
Z = Z / sum(Z(:));
 
Z = Z / max(Z(:));  
cauchyImage = uint8(Z * 255);
Z=Z/(sum(sum(Z)));
filename = 'dualcauchy.mat';
save(filename, 'Z');
figure('Position', [100, 100, 800, 800]);
imshow(cauchyImage);
%title('Test Measure: 2D Cauchy Distribution (1024x1024)', 'FontSize', 14);
colormap('parula');  
colorbar;
% imwrite(cauchyImage, 'cauchy_1024.png');