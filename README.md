# Processamento-Digital-de-Imagens

## Lista 2 Códigos Matlab

* Ex1:
~~~Matlab
%limpa tela e fecha a janela
close all;
clear;
clc;

%Declaração Sinal  
fs = 200; f = 10;
time = 0 : 1/fs : (1 - 1/fs);

%Plotando o sinal Original
k = 0;
x = sin(2 * pi * f * time) + k * randn(1, fs);
subplot(2,2,1);
plot(x)
title('Sinal X')
ylabel('Valor')
xlabel('Frequência angular')

%Plotando o FFT para k = 0.5
k = 0.5
y = fft(x);
m = abs(y); 
subplot(2,2,2);
plot(m)
title('FFT do Sinal X para k = 0.5')
ylabel('Valor')
xlabel('Frequência angular')

%Plotando o FFT para k = 1.5
k = 1.5;
x = sin(2 * pi * f * time) + k * randn(1, fs);
yy = fft(x);
m = abs(yy); 
subplot(2,2,3);
plot(m)
title('FFT do Sinal X para k = 1.5')
ylabel('Valor')
xlabel('Frequência angular')

%Plotando o FFT para k = 3
k = 3;
x = sin(2 * pi * f * time) + k * randn(1, fs);
y = fft(x);
m = abs(y); 
subplot(2,2,4);
plot(m)
title('FFT do Sinal X para k = 3')
ylabel('Valor')
xlabel('Frequência angular')

~~~
* Ex2:
~~~Matlab
%limpa tela e fecha a janela
close all;
clear;
clc;

Npoints   = 251; % BE CAREFUL WITH THIS NUMBER!
                   % DTFT is O(n^2)
                   
%preparação do sinal
t = 0 : 0.001 : 0.25;
Signal = sin(2 * pi * 50 * t) + sin(2 * pi * 120 * t)

% Começo do timer pro DTFT
tic;

% DTFT calculo força bruta
FourierTime = (0:Npoints-1)/Npoints;
FourierCoefs   = zeros(size(Signal));

for i = 1:Npoints
    ComplexSineWave = exp( -1j*2*pi*(i-1)*FourierTime );
    FourierCoefs(i) = sum( Signal.*ComplexSineWave );
end

% Timer pro DTFT
z(1) = toc;

% começo do timer pro FFT
tic;
FourierCoefsFFT = fft(Signal);
z(2) = toc;

%Fazendos o plot dos dados para comparação futura
subplot(2,1,1)
a = abs(FourierCoefs)
plot(a)
title('DFT')
xlabel('tempo')
ylabel('valor')

subplot(2,1,2)
b = abs(FourierCoefsFFT)
plot(b)
title('FFT')
xlabel('tempo')
ylabel('valor')
~~~

* Ex3:
~~~Matlab
%limpa tela e fecha a janela
close all;
clear;
clc;

%Declaração Sinal  
y(1 : 10) = 1;
y(248 : 256) = 1;


%Plotando o sinal Original
subplot(2,1,1);
plot(y)
title('Sinal Y')
ylabel('Valor')
xlabel('Tempo')

%Plotando o sinal IFFT
subplot(2,1,2);
a = ifft(y)
b = abs(a)
plot(b)
title('Transformada Inversa de Fourier do Sinal Y')
ylabel('Valor')
xlabel('Tempo')


~~~
