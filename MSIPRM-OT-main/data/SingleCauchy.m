%clear; close all; clc;
 
imageSize = 1024;  
k=imageSize/1024;
mu = [imageSize/2, imageSize/2]+[100,100]*k;  
gamma = [600, 1000]*k;               
rho = 0.5;

% imageSize = 1024;  
% k=imageSize/1024;
% mu = [imageSize/2, imageSize/2];  
% gamma = [300, 1000]*k;               
% rho = 0.5;                        


[X, Y] = meshgrid(1:imageSize, 1:imageSize);

x = (X - mu(1)) / gamma(1);
y = (Y - mu(2)) / gamma(2);

Z = 1 ./ (2*pi*sqrt(1-rho^2)) .* ...
    1 ./ (1 + (x.^2 - 2*rho*x.*y + y.^2)/(1-rho^2)).^(3/2);

Z = Z / max(Z(:));  
cauchyImage = uint8(Z * 255);
Z=Z/(sum(sum(Z)));
filename = 'Singlecauchy.mat';
save(filename, 'Z');
figure('Position', [100, 100, 800, 800]);
imshow(cauchyImage);
%title('Test Measure: 2D Cauchy Distribution (1024x1024)', 'FontSize', 14);
colormap('parula');  
colorbar;
% imwrite(cauchyImage, 'cauchy_1024.png');