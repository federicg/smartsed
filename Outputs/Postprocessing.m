close all
clear
clc

folder="0";


%% Display constants

filename = char( folder+"/DEM.asc");
[b, R] = arcgridread(filename);

filename = char( folder+"/basin_mask.asc");
[basin_mask, R] = arcgridread(filename);

filename = char( folder+"/soilMoistureRetention.asc");
[soilMoistureRetention, R] = arcgridread(filename);

filename = char( folder+"/slope_x.asc");
[slope_x, R] = arcgridread(filename);

filename = char( folder+"/slope_y.asc");
[slope_y, R] = arcgridread(filename);

filename = char( folder+"/d_10.asc");
[d_10, R] = arcgridread(filename);

filename = char( folder+"/d_90.asc");
[d_90, R] = arcgridread(filename);

filename = char( folder+"/k_c.asc");
[k_c, R] = arcgridread(filename);

filename = char( folder+"/excluded_ids_bool.asc");
[exc_ids_bool, R] = arcgridread(filename);

filename = char( folder+"/excluded_ids_pour.asc");
[exc_ids_pour, R] = arcgridread(filename);


S_x = ( slope_x(:,1:end-1) + slope_x(:,2:end) ) * .5;
S_y = ( slope_y(1:end-1,:) + slope_y(2:end,:) ) * .5;
S = (S_x.^2 + S_y.^2).^0.5;


cellsize = R(2,1);

fprintf("cellsize: ")
disp(cellsize)

%%%%% build coordinate matrices

Nx = size(b,2);
Ny = size(b,1);
j = 1:Nx;
i = 1:Ny;
[jj,ii] = meshgrid(j,i);
XX = (jj-1).*cellsize;
YY = - (ii-1) .* cellsize + (Ny-1).*cellsize;

disp('# excluded cells, ')
disp(length(exc_ids_bool(exc_ids_bool==1)))

disp('% excluded cells over basin cells, ')
disp(length(exc_ids_bool(exc_ids_bool==1))/numel(basin_mask(basin_mask==1))*100)

basin_mask(~basin_mask) = NaN;
exc_ids_bool(exc_ids_bool==1)=NaN;
exc_ids_bool(exc_ids_bool==0)=1;

exc_ids_pour(exc_ids_pour==-1)=NaN;
exc_ids_pour_i = floor(exc_ids_pour./Nx); 
exc_ids_pour_j = mod(exc_ids_pour, Nx);

ii=ii-1;
jj=jj-1;
II = exc_ids_pour_i-ii;
JJ = exc_ids_pour_j-jj;


figure()
colormap( copper )
mesh(XX/1000,YY/1000,b.*basin_mask)
xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
colorbar

figure()
colormap( copper )
mesh(XX/1000,YY/1000,exc_ids_bool.*b.*basin_mask)
xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
colorbar

figure()
colormap( copper )
mesh(XX/1000,YY/1000,exc_ids_bool.*b.*basin_mask)
xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
colorbar


figure()
colormap( copper )
imagesc(exc_ids_bool)
hold on
quiver(jj+1,ii+1,JJ,II,0);
xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
colorbar

figure()
colormap( copper )
contourf(S.*basin_mask)
axis ij
hold on
quiver(jj+1,ii+1,JJ,II,0);
xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
colorbar


figure()
colormap( copper )
contourf(XX/1000,YY/1000,soilMoistureRetention.*basin_mask)
xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
title('max soil mois. ret.')
colorbar

% figure()
% colormap( copper )
% contourf(slope_x)
% xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
% ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
% title('slope_x')
% colorbar
% 
% figure()
% colormap( copper )
% contourf(slope_y)
% xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
% ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
% title('slope_y')
% colorbar
%
figure()
% colormap( copper )
contourf(XX/1000,YY/1000,exc_ids_bool.*S.*basin_mask)
xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
title('slope: $\|\nabla b\|$','fontsize',10,'interpreter','latex')
colorbar
%
figure()
colormap( copper )
contourf(XX/1000,YY/1000,d_10.*basin_mask)
xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
title('d_{10}')
colorbar

figure()
colormap( copper )
contourf(XX/1000,YY/1000,d_90.*basin_mask)
xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
title('d_{90}')
colorbar

figure()
colormap( winter )
contourf(XX/1000,YY/1000,k_c.*basin_mask)
xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
title('k_c')
colorbar



%% Display time variable solutions

numberSimDays=input('Insert below the number of simulated days \n');



%%

filename = char( folder+"/u_" + string(0) + ".asc" );
[u, R] = arcgridread(filename);

filename = char( folder+"/v_" + string(0) + ".asc");
[v, R] = arcgridread(filename);

% center cell velocities
u = ( u(:,1:end-1) + u(:,2:end) ) * .5;
v = ( v(1:end-1,:) + v(2:end,:) ) * .5;
vel(:,:,1) = (u.^2 + v.^2).^0.5;
% vel=zeros(size(b,1),size(b,2),numberSimDays);
video = VideoWriter('vel.avi');
video.FrameRate = 10; %set your frame rate
open(video);
figure()
contourf(XX/1000,YY/1000,vel(:,:,1).*basin_mask)
xlabel('$x \; (km)$','fontsize',10,'interpreter','latex')
ylabel('$y \; (km)$','fontsize',10,'interpreter','latex')
colormap(winter)
lim = caxis;
colorbar
set(gca,'nextplot','replacechildren');
for n = 1:1:numberSimDays
    n
    filename = char( folder+"/u_" + string(n) + ".asc" );
    [u, R] = arcgridread(filename);
    
    filename = char( folder+"/v_" + string(n) + ".asc");
    [v, R] = arcgridread(filename);
    
    % center cell velocities
    u = ( u(:,1:end-1) + u(:,2:end) ) * .5;
    v = ( v(1:end-1,:) + v(2:end,:) ) * .5;
    vel(:,:,n) = (u.^2 + v.^2).^0.5;
    contourf(XX/1000,YY/1000,vel(:,:,n).*basin_mask)%,linspace(lim(1),lim(2),1000))
    %zlim([lim(1), lim(2)]);
    drawnow;
    frame = getframe;
    writeVideo(video,frame);
    
    
end
close(video);

% 
%%
% 
video = VideoWriter('vel_vect.avi');
video.FrameRate = 10; %set your frame rate
open(video);
figure()
[M,c]=contour(XX/1000,YY/1000,b.*basin_mask);
c.LineWidth = 3;
hold on
quiver(XX/1000,YY/1000,u(1:1:end,1:1:end)*0, ...
       -v(1:1:end,1:1:end)*0,'k','AutoScaleFactor',3);
xlabel('$x \; (km)$','fontsize',10,'interpreter','latex')
ylabel('$y \; (km)$','fontsize',10,'interpreter','latex')
colormap(copper)
colorbar
set(gca,'nextplot','replacechildren');
for n = 1:numberSimDays
    n
    filename = char( folder+"/u_" + string(n) + ".asc" );
    [u, R] = arcgridread(filename);
    
    filename = char( folder+"/v_" + string(n) + ".asc");
    [v, R] = arcgridread(filename);
    
    % center cell velocities
    u = ( u(:,1:end-1) + u(:,2:end) ) * .5.*basin_mask;
    v = ( v(1:end-1,:) + v(2:end,:) ) * .5.*basin_mask;
    
    [M,c]=contour(XX/1000,YY/1000,b.*basin_mask);
    colorbar
    c.LineWidth = 3;
    hold on
    quiver(XX/1000,YY/1000,u(1:1:end,1:1:end), ...
        -v(1:1:end,1:1:end),'k','AutoScaleFactor',3);
    xlabel('$x \; (km)$','fontsize',10,'interpreter','latex')
    ylabel('$y \; (km)$','fontsize',10,'interpreter','latex')
    colormap(copper)
    hold off
%     pause
    frame = getframe;
    
    writeVideo(video,frame);
%     
    
end
close(video);



%%

filename = char( folder+"/H_" + string(0) + ".asc");
[H(:,:,1), R] = arcgridread(filename);
video = VideoWriter('H.avi');
video.FrameRate = 10; %set your frame rate
open(video);
figure()
contourf(XX/1000,YY/1000,H(:,:,1))
xlabel('$x \; (km)$','fontsize',10,'interpreter','latex')
ylabel('$y \; (km)$','fontsize',10,'interpreter','latex')
colorbar
colormap(winter)
% pause
lim = caxis;
set(gca,'nextplot','replacechildren');
for n = 1:numberSimDays
    filename = char( folder+"/H_" + string(n) + ".asc");
    [H(:,:,n), R] = arcgridread(filename);
    n
    contourf(XX/1000,YY/1000,H(:,:,n).*basin_mask)%,linspace(lim(1),lim(2),30))
%     return
    %     caxis(lim)
%     colorbar
% pause
%     zlim([lim(1), lim(2)]);
    drawnow;
    frame = getframe;
    
    writeVideo(video,frame);
    %pause
end
close(video);

figure()
colormap( copper )
mesh(XX/1000,YY/1000,b.*basin_mask)
hold on
surf(XX/1000,YY/1000,(b+H(:,:,end)).*basin_mask)
xlabel('$x\;(km)$','fontsize',10,'interpreter','latex')
ylabel('$y\;(km)$','fontsize',10,'interpreter','latex')
colorbar


%%

w_cum=zeros(size(b,1),size(b,2),numberSimDays);
video = VideoWriter('w_cum.avi');
video.FrameRate = 1; %set your frame rate
open(video);
contourf(XX/1000,YY/1000,w_cum(:,:,1))
xlabel('$x \; (km)$','fontsize',10,'interpreter','latex')
ylabel('$y \; (km)$','fontsize',10,'interpreter','latex')
colorbar
set(gca,'nextplot','replacechildren');
for n = 1:numberSimDays
    filename = char( folder+"/w_cum_" + string(n) +".asc");
    [w_cum(:,:,n), R] = arcgridread(filename);
    n
    
    
    contourf(XX/1000,YY/1000,w_cum(:,:,n).*basin_mask)
    colorbar
    frame = getframe;
    writeVideo(video,frame);
end
close(video);



%%

h_G=zeros(size(b,1),size(b,2),numberSimDays);
video = VideoWriter('h_G.avi');
video.FrameRate = 1; %set your frame rate
open(video);
contourf(XX/1000,YY/1000,h_G(:,:,1))
% imagesc(h_G(:,:,1).*basin_mask)
xlabel('$x \; (km)$','fontsize',10,'interpreter','latex')
ylabel('$y \; (km)$','fontsize',10,'interpreter','latex')
colorbar
set(gca,'nextplot','replacechildren');
for n = 1:numberSimDays
    filename = char( folder+"/hG_" + string(n) +".asc");
    [h_G(:,:,n), R] = arcgridread(filename);
    
    
    
    contourf(XX/1000,YY/1000,h_G(:,:,n).*basin_mask)
    %imagesc(h_G(:,:,n).*basin_mask)
    %colorbar
    frame = getframe;
    writeVideo(video,frame);
end
close(video);

%%

ET=zeros(size(b,1),size(b,2),numberSimDays);
video = VideoWriter('ET.avi');
video.FrameRate = 1; %set your frame rate
open(video);
contourf(XX/1000,YY/1000,ET(:,:,1))
% imagesc(h_G(:,:,1).*basin_mask)
xlabel('$x \; (km)$','fontsize',10,'interpreter','latex')
ylabel('$y \; (km)$','fontsize',10,'interpreter','latex')
colorbar
set(gca,'nextplot','replacechildren');
for n = 1:numberSimDays
    filename = char( folder+"/ET_" + string(n) +".asc");
    [ET(:,:,n), R] = arcgridread(filename);
    
    
    
    contourf(XX/1000,YY/1000,ET(:,:,n).*basin_mask)
    %imagesc(h_G(:,:,n).*basin_mask)
    %colorbar
    frame = getframe;
    writeVideo(video,frame);
end
close(video);

%%

h_sd=zeros(size(b,1),size(b,2),numberSimDays);
video = VideoWriter('h_sd.avi');
video.FrameRate = 10; %set your frame rate
open(video);
contour(XX/1000,YY/1000,h_sd(:,:,1))
xlabel('$x \; (km)$','fontsize',10,'interpreter','latex')
ylabel('$y \; (km)$','fontsize',10,'interpreter','latex')
colorbar
colormap( copper )
set(gca,'nextplot','replacechildren');
for n = 1:1:numberSimDays
    n
    filename = char( folder+"/hsd_" + string(n) +".asc");
    [h_sd(:,:,n+1), R] = arcgridread(filename);
    [M,c]=contour(XX/1000,YY/1000,b.*basin_mask);
    colorbar
    c.LineWidth = 3;
    hold on
    
    
    contour(XX/1000,YY/1000,h_sd(:,:,n).*basin_mask)
    hold off
    frame = getframe;
    writeVideo(video,frame);
end
close(video);

%%

h_sn=zeros(size(b,1),size(b,2),numberSimDays);
video = VideoWriter('h_sn.avi');
video.FrameRate = 1; %set your frame rate
open(video);
contourf(XX/1000,YY/1000,h_sn(:,:,1))
xlabel('$x \; (km)$','fontsize',10,'interpreter','latex')
ylabel('$y \; (km)$','fontsize',10,'interpreter','latex')
colorbar
set(gca,'nextplot','replacechildren');
for n = 0:numberSimDays
    filename = char( folder+"/hsn_" + string(n) +".asc");
    [h_sn(:,:,n+1), R] = arcgridread(filename);
    
    
    
    contourf(XX/1000,YY/1000,h_sn(:,:,n+1).*basin_mask)
    colorbar
    frame = getframe;
    writeVideo(video,frame);
end
close(video);



%%

numberPrec=input('Insert below the number \n');

f=zeros(size(b,1),size(b,2),numberPrec);
p=zeros(size(b,1),size(b,2),numberPrec);


video = VideoWriter('f.avi');
video.FrameRate = 1; %set your frame rate
open(video);
contourf(XX/1000,YY/1000,f(:,:,1))
xlabel('$x \; (km)$','fontsize',10,'interpreter','latex')
ylabel('$y \; (km)$','fontsize',10,'interpreter','latex')
colorbar
set(gca,'nextplot','replacechildren');
for n = 1:numberPrec
    filename = char( folder+"/f_" + string(n) +".asc");
    [f(:,:,n+1), R] = arcgridread(filename);
    
    
    contourf(XX/1000,YY/1000,f(:,:,n+1).*basin_mask)
    %colorbar
    frame = getframe;
    writeVideo(video,frame);
end
close(video);


%%
video = VideoWriter('p.avi');
video.FrameRate = 1; %set your frame rate
open(video);
contourf(XX/1000,YY/1000,p(:,:,1))
xlabel('$x \; (km)$','fontsize',10,'interpreter','latex')
ylabel('$y \; (km)$','fontsize',10,'interpreter','latex')
colorbar
set(gca,'nextplot','replacechildren');
for n = 1:numberPrec
    filename = char( folder+"/p_" + string(n) +".asc");
    [p(:,:,n+1), R] = arcgridread(filename);
    
    
    contourf(XX/1000,YY/1000,p(:,:,n).*basin_mask)
    colorbar
    frame = getframe;
    writeVideo(video,frame);
end
close(video);



%% Setup the Import Options and import the data


opts = spreadsheetImportOptions("NumVariables", 4);

% Specify sheet and range
opts.Sheet = "Evento Novembre 2018";
opts.DataRange = "B554:E1595";

% Specify column names and types
opts.VariableNames = ["Ora", "ViaCarloPortaLIVELLOIDROmediadelperiodo1hcm", "ViaCarloPortaLIVELLOIDROmassimodelperiodo1hcm", "ViaCarloPortaPIOGGIAtotale1hmm"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Import the data
datipioggiaportataS1 = readtable("caldone_carloPorta/dati_pioggia_portata.xlsx", opts, "UseExcel", false);

% Convert to output type
datipioggiaportataS1 = table2array(datipioggiaportataS1);

% Clear temporary variables
clear opts


% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 1);

% Specify sheet and range
opts.Sheet = "Evento Novembre 2018";
opts.DataRange = "A554:A1595";

% Specify column names and types
opts.VariableNames = "Data";
opts.VariableTypes = "datetime";

% Specify variable properties
opts = setvaropts(opts, "Data", "InputFormat", "");

% Import the data
times = readtable("../../dati_pioggia_portata.xlsx", opts, "UseExcel", false);

% Clear temporary variables
clear opts


times=times(1:end,"Data").Variables;
times.Hour=floor(datipioggiaportataS1(:,1)*24);
times.Minute=round((datipioggiaportataS1(:,1)*24-floor(datipioggiaportataS1(:,1)*24))*60);


waterHeight=abs(datipioggiaportataS1(:,2)); % c'è un valore acqua negativo
waterHeight=waterHeight/100; % convert to meters

figure()
plot(times, waterHeight)
hold on


waterHeight_num  = load('0/waterSurfaceHeight_1.txt');
times_sim = waterHeight_num(:,2)/3600/24+times(1)-4;

plot(times_sim,waterHeight_num(:,1),'o')
xlim([times(1),times(end)])
xlabel('$time\:(sec.)$','interpreter','latex')
ylabel('$H\:(m)$','interpreter','latex')



%% portata via Carlo Porta

Q = @(H) 25.578.*H.^(3.098);

folder=["45218.hpc.mate.polimi.it/Outputs/0", "45216.hpc.mate.polimi.it/Outputs/0", "45217.hpc.mate.polimi.it/Outputs/0"]; %["44765.hpc.mate.polimi.it/Outputs/0", "44771.hpc.mate.polimi.it/Outputs/0"];
cellsize_=[20,35,50];

figure()
plot(times, Q(waterHeight), '--b')
hold on
for i = 1:length(folder)
    waterHeight_num  = readmatrix("../"+folder(i)+"/waterSurfaceMassFlux.txt");
    i
    dt_numm=waterHeight_num(1);
    waterHeight_num=waterHeight_num(2:end);
    
    days_simulated_numm=dt_numm*length(waterHeight_num)/(3600*24);
    
    times_simulation_num =linspace(times(1),times(1)+days_simulated_numm,length(waterHeight_num))';
    
    waterHeight_num(logical((times_simulation_num.Month>10) ...
        .* (times_simulation_num.Day>18))) = nan;
    
    plot(times_simulation_num, waterHeight_num*cellsize_(i),'o')
    hold on
end
hold off
legend('real','20 (m)', '35 (m)','50 (m)')
xlabel('$time\:(sec.)$','interpreter','latex')
ylabel('$Q\:(m^3/sec.)$','interpreter','latex')

figure()
plot(times_simulation_num, waterHeight_num*cellsize,'o')


%% read matrix and vector
clc


load ../0/A_1.txt
A_1 = spconvert(A_1);
%issymmetric(A)
try chol(A_1);
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end


load ../0/b_1.txt
load ../0/H_1.txt



figure()
plot(H_1-A_1\b_1)

figure()
plot(H_1)

figure()
plot(A_1\b_1)
title('matlab')



load ../0/A_2.txt
A_2 = spconvert(A_2);
%issymmetric(A)
try chol(A_2);
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end


load ../0/b_2.txt
load ../0/H_2.txt



figure()
plot(H_2-A_2\b_2)

figure()
plot(H_2)

figure()
plot(A_2\b_2)
title('matlab')


% b_2=H_1
figure()
plot(b_2-H_1)




% %%
% 
% R = @(ep1,ep2) (sqrt(ep1)-sqrt(ep2)) ./ (sqrt(ep1)+sqrt(ep2));
% 
% ep1=9;
% ep2=linspace(1,80,100);
% R_ep2 = R(ep1,ep2);
% 
% figure()
% plot(ep2,R_ep2)



%%

load('0/waterSurfaceHeight.txt');

H = waterSurfaceHeight(:,1);
time = waterSurfaceHeight(:,2);


start_day = datetime(2018,10,22);
t1 = start_day + time/(3600*24);

% figure()
% plot(time, H)
% return


opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "Evento Novembre 2018";
opts.DataRange = "A2:E1595";

% Specify column names and types
opts.VariableNames = ["Data", "Ora", "ViaCarloPortaLIVELLOIDROmediaDelPeriodo1hcm", "ViaCarloPortaLIVELLOIDROmassimoDelPeriodo1hcm", "ViaCarloPortaPIOGGIAtotale1hmm"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Data", "InputFormat", "");

% Import the data
datipioggiaportata = readtable("/Users/federicog/Documents/smartsed/github_repo/smartsed/Outputs/57538.hpc.mate.polimi.it/Outputs/caldone_carloPorta/dati_pioggia_portata.xlsx", opts, "UseExcel", false);
clear opts




t2 = table2array(datipioggiaportata(:,1)) + table2array(datipioggiaportata(:,2));
hh = table2array(datipioggiaportata(:,3))/100;


figure()
plot(t1,H,'o--') 
hold on
plot(t2, hh,'o--') 
xlim([datetime(2018,10,26),t2(end)])
ylabel('$H\:(m)$','interpreter','latex')



% portata
Q = @(H) 25.578.*H.^(3.098);


load('0/waterSurfaceMassFlux.txt');
QQ = waterSurfaceHeight(:,1)*50;

figure()
plot(t1,Q(H),'o--')
hold on
plot(t2,Q(hh),'o--')
xlim([datetime(2018,10,26),t2(end)])
ylabel('$Q\:(m^3/sec.)$','interpreter','latex')


