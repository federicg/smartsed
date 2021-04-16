close all
clear
clc



%% Preprocessing of Smart-sed

[clay, R] = readgeoraster('clay.asc');
[sand, R] = readgeoraster('sand.asc');
[clc,  R] = readgeoraster('CLC.asc');

Nx = R.RasterSize(2);
Ny = R.RasterSize(1);

dx = R.CellExtentInWorldX;


file_id   = fopen('clay_id.asc', 'w');

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
        fprintf(file_id, '%12.8f ', 0.3);
    end
    fprintf(file_id, '\n');
end
fclose(file_id);



file_id   = fopen('sand_id.asc', 'w');

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
        fprintf(file_id, '%12.8f ', 0.3);
    end
    fprintf(file_id, '\n');
end
fclose(file_id);




file_id   = fopen('CLC_id.asc', 'w');

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
        fprintf(file_id, '%12.8f ', 141);
    end
    fprintf(file_id, '\n');
end
fclose(file_id);




