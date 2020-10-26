function [G,f] = dft(g,t,Nf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRASFORMATA DISCRETA DI FOURIER 
% g = segnale da trasformare
% t = asse dei tempi
% Nf = numero dei punti in frequenza da calcolare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The execution time for fft depends on the length of the transform. 
% The time is fastest for powers of two and almost as fast for lengths 
% that have only small prime factors. The time is typically several 
% times slower for lengths that are prime, or which have large prime factors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(exist('Nf'))
    Nf = length(t);
end

g = g(:); % sempre interpretato come vettore colonna

G = fft(g,Nf);
G = fftshift(G);  % per centrare attorno a zero

% asse delle frequenze
dt = t(2)-t(1);
if mod(Nf,2) % numero dispari di campioni
    f = (-(Nf-1)/2:(Nf-1)/2)/Nf/dt;
else   % numero pari di campioni
    f = (-Nf/2:Nf/2-1)/Nf/dt;
end
f = f';

% compensazione del ritardo della fft
fase = -2*pi*f(:)*t(1);

G = G.*exp(1i*fase(:));

end