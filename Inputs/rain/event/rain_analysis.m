%close all
clear
clc


%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 3);

% Specify sheet and range
opts.Sheet = "Anno 2018";
opts.DataRange = "A2:C8818";

% Specify column names and types
opts.VariableNames = ["Data", "Ora", "ViaCarloPortaPIOGGIAtotale1hmm"];
opts.VariableTypes = ["double", "double", "double"];

% Import the data
pioggiaS1 = readtable("/Users/federicog/Documents/my_smartsed/ReleaseLike/Inputs/rain/event/pioggia.xlsx", opts, "UseExcel", false);

% Convert to output type
pioggiaS1 = table2array(pioggiaS1);

% Clear temporary variables
clear opts


%% Find ww s.t. is maximum among feasibles and corrcoeff sia >= threshold

ww=100;

A=pioggiaS1(:,3);
B=A;
B(1:ww)=0;

dt=1/24; % giorni
t=(1:length(A))*dt;
[GA,f]=dft(A-mean(A),t);
[GB,f]=dft(B-mean(B),t);


figure()
stem(f,abs(GA))

figure()
stem(f,abs(GB))

figure()
stem(f,angle(GA))

figure()
stem(f,angle(GB))


figure()
scatter(abs(GB),abs(GA),'b')
hold on
plot(abs(GB),abs(GB),'--r*')

figure()
scatter(angle(GB),angle(GA),'b')
hold on
plot(angle(GB),angle(GB),'--r*')

