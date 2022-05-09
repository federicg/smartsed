close all
clear
clc



%% Preprocessing of Smart-sed



% In these file an idealized orography is set and output as ASCII files


% set the length of the domain in meters
dx = 5;
dy = dx;


Nx = 2000;
Ny = 2000; %1829


Lx = Nx*dx;
Ly = Ny*dy;




%% Set initial condition on H


radius = 0.25*Lx; %[m]
j = 1:Nx;
i = 1:Ny;
[jj,ii] = meshgrid(j,i);
XX = (jj-1).*dx;
YY = - (ii-1) .* dy + (Ny-1).*dy;
r = sqrt(( YY - Ly*0.5 ).^2 + ( XX - Lx*0.5 ).^2 );
% H = 1.e-3.*ones(size(XX));
% H = (0.02.*(cos(r * pi * 0.5 ./ radius)).^3).*(r<=radius);
H = 1 + 1*exp(-0.5*( (XX-Lx/2).^2+(YY-Ly/2).^2 )/(0.2*Lx)^2);
% H = 0.02.*(r<=radius);
% H = 0.02.*(XX<=Lx/2);
oro = (2 + Ly*.0005+ .0005 .* XX - .0005 .* YY);
%oro = 15-80*((XX-Lx*0.5).^2 +(YY-Ly*0.5).^2)./Lx^2 +300*((XX-Lx*0.5).^4 +(YY-Ly*0.5).^4)./Lx^4;
%oro = oro*0;


if ( sum(sum(H < 0)) ~= 0 )
    disp("one H is negative!!");
    %return;
end





%%
figure()
colormap(winter)
mesh(XX/1000,YY/1000,H)
xlabel('south')
ylabel('west')
colorbar
axis equal
title('Initial surface water layer depth [m]')


figure()
colormap(winter)
contourf(XX/1000,YY/1000,oro)
xlabel('south')
ylabel('west')
colorbar
axis equal



% % for plot unless use of imagesc --> H_
% for i = 1:Ny
%     H_(i,:) = H(end-i+1,:);
%     b_(i,:) = b(end-i+1,:);    
% end
% 
% 
% 
% figure(3)
% surf(H_)
% colormap(jet)
% xlabel('south')
% ylabel('west')


%% Print output files



file_id   = fopen('H.asc', 'w');


%
fprintf(file_id, '%5s ', 'ncols');
fprintf(file_id, '%u\n', Nx);

%
fprintf(file_id, '%5s ', 'nrows');
fprintf(file_id, '%u\n', Ny);

%
fprintf(file_id, '%9s ', 'xllcorner');
fprintf(file_id, '%12.8f\n', 528669.646615625359);

%
fprintf(file_id, '%9s ', 'yllcorner');
fprintf(file_id, '%12.8f\n', 5076881.926118221134);

%
fprintf(file_id, '%8s ', 'cellsize');
fprintf(file_id, '%12.8f\n', dx);

%
fprintf(file_id, '%12s ', 'NODATA_value');
fprintf(file_id, '%12.8f\n', -9999);




% outputstr = ['%' num2str(size(H,1)) 'f '];
% outputstr = repmat(outputstr, 1, size(H,2));
% outputstr = [outputstr '\n'];
% fprintf(file_id, outputstr, H.');
% fclose(file_id);

for i = 1:Ny
    for j = 1:Nx
        fprintf(file_id, '%12.8f ', H(i,j));
    end
    fprintf(file_id, '\n');
end
fclose(file_id);


% %%
% matrix = magic(4) % example matrix
% [mrows, ncols] = size(matrix)
% outputstr = ['%' num2str(mrows) 'f ']
% outputstr=repmat(outputstr, 1, ncols)
% outputstr = [outputstr '\n']
% fprintf(file_id,outputstr, matrix.')
%% Print output files



file_id   = fopen('DEMIdeal.asc', 'w');


%
fprintf(file_id, '%5s ', 'ncols');
fprintf(file_id, '%u\n', Nx);

%
fprintf(file_id, '%5s ', 'nrows');
fprintf(file_id, '%u\n', Ny);

%
fprintf(file_id, '%9s ', 'xllcorner');
fprintf(file_id, '%12.8f\n', 528669.646615625359);

%
fprintf(file_id, '%9s ', 'yllcorner');
fprintf(file_id, '%12.8f\n', 5076881.926118221134);

%
fprintf(file_id, '%8s ', 'cellsize');
fprintf(file_id, '%12.8f\n', dx);

%
fprintf(file_id, '%12s ', 'NODATA_value');
fprintf(file_id, '%12.8f\n', -9999);





for i = 1:Ny
    for j = 1:Nx
        fprintf(file_id, '%12.8f ', oro(i,j));
    end
    fprintf(file_id, '\n');
end
fclose(file_id);

    

%%
file_id   = fopen('mask.asc', 'w');


%
fprintf(file_id, '%5s ', 'ncols');
fprintf(file_id, '%u\n', Nx);

%
fprintf(file_id, '%5s ', 'nrows');
fprintf(file_id, '%u\n', Ny);

%
fprintf(file_id, '%9s ', 'xllcorner');
fprintf(file_id, '%12.8f\n', 528669.646615625359);

%
fprintf(file_id, '%9s ', 'yllcorner');
fprintf(file_id, '%12.8f\n', 5076881.926118221134);

%
fprintf(file_id, '%8s ', 'cellsize');
fprintf(file_id, '%12.8f\n', dx);

%
fprintf(file_id, '%12s ', 'NODATA_value');
fprintf(file_id, '%12.8f\n', -9999);




a = 0;
for i = 1:Ny
    for j = 1:Nx
        if ( i > a && i < Ny-a+1 && j > a && j < Nx-a+1 )
            fprintf(file_id, '%12.8f ', 1);
        else
            fprintf(file_id, '%12.8f ', 1);
        end
    end
    fprintf(file_id, '\n');
end
fclose(file_id);
return
%%
basin = geotiffread('Mask_bin.tif');

figure()
contourf(XX/1000,YY/1000,basin)
colorbar

% fun=YY/1000
% quad=basin.*(fun>-1.5);
fun=-4+XX/1000+YY/1000;
quad=basin.*(fun>0);

figure()
contourf(XX/1000,YY/1000,quad)
colorbar



file_id   = fopen('Mask_bin_my.asc', 'w');


%
fprintf(file_id, '%5s ', 'ncols');
fprintf(file_id, '%u\n', Nx);

%
fprintf(file_id, '%5s ', 'nrows');
fprintf(file_id, '%u\n', Ny);

%
fprintf(file_id, '%9s ', 'xllcorner');
fprintf(file_id, '%12.8f\n', 528669.646615625359);

%
fprintf(file_id, '%9s ', 'yllcorner');
fprintf(file_id, '%12.8f\n', 5076881.926118221134);

%
fprintf(file_id, '%8s ', 'cellsize');
fprintf(file_id, '%12.8f\n', dx);

%
fprintf(file_id, '%12s ', 'NODATA_value');
fprintf(file_id, '%12.8f\n', -9999);



for i = 1:Ny
    for j = 1:Nx
        fprintf(file_id, '%12.8f ', quad(i,j));
    end
    fprintf(file_id, '\n');
end
fclose(file_id);



