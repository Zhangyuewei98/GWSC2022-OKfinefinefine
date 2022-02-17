clc;clear;close all;

addpath CODES\
ffparams = struct('rmin',-5,...
                     'rmax',5 ...
                  );
[X,Y] = meshgrid(linspace(0, 1, 500),linspace(0, 1, 500));
Z = crcbpsotestfunc([X(:)'; Y(:)']', ffparams);
Z = reshape(Z, 500, 500);
Z_min = min(min(Z));
Z_max = max(max(Z));

figure
surf(X, Y, Z, 'LineStyle', 'none');
axis xy;
xlabel('x_1');
ylabel('x_2');
zlabel('f(x_1, x_2)');
saveas(gcf,'test_surf','jpg');