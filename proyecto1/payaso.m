% Equipo: Mariana, Rom�n, Santiago
% 8 de Octubre del 2020
% Proyecto 1 Optimizaci�n
% -------------------------------------------------------------------------
% Obtenci�n de la descomposici�n de im�genes en dos matrices W, H
% v�a el m�todo de descenso en dos pasos a partir del m�todo de puntos
% interiores.
% -------------------------------------------------------------------------
% Configuraci�n inicial
clear; clc; close all; warning('off');

% Obtener imagen en el ambiente local
load('clown')
colormap('gray')

% obtenci�n de W, H para k = 5, 20, 30, 60, 80
k = [5; 20; 30; 60; 80];

% metodo de optimizaci�n quadprog
% salvar estructuras
WH_qp = zeros(200, 320, 5);
time_qp = zeros(5);
normasQP = zeros(5,1);

% metodo
for i = 1:numel(k)
    t0 = cputime;
    [W, H] = descenso2pasos_qp(X, k(i));
    WH_qp(:,:,i) = W*H;
    time_cpu = cputime - t0;
    time_qp(i) = time_cpu;
    normasQP(i) = norm(X-W*H, 'fro');
end 

WH_Nuestro = zeros(200, 320, 5);
time_Nuestro = zeros(5, 1);
normas = zeros(5, 1);
for i = 1:numel(k)
    t0 = cputime;
    [W, H] = descenso2pasos(X, k(i));
    WH_Nuestro(:,:,i) = W*H;
    time_cpu = cputime - t0;
    time_Nuestro(i) = time_cpu;
    normas(i) = norm(X-W*H, 'fro');
end

% Hacemos las comparaciones

figure(1)
scatter(k, time_Nuestro, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'b')
hold on
scatter(k, time_qp(:,1), 'MarkerEdgeColor', [0 0 0.8], 'MarkerFaceColor', [0.6 0.2 0.9])
legend('Nuestro metodo', 'QuadProg')
title('Tiempo que tarda nuestro programa vs tamaño de k')
xlabel('k')
ylabel('segundos')

figure(2)
scatter(k, normas, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'b')
hold on
scatter(k, normasQP, 'MarkerEdgeColor', [0 0 0.8], 'MarkerFaceColor', [0.6 0.2 0.9])
legend('Nuestro metodo', 'QuadProg')
title('Norma Frobenius de X-W*H vs tamaño de k')
xlabel('k')
ylabel('1.0e-05 * Norma Frobenius de X-W*H')
ylim([0 3.5e+03])


% Y sacamos las imagenes

%Quadprog
figure(3)
colormap('gray')
image(WH_qp(:,:,1))

figure(4)
colormap('gray')
image(WH_qp(:,:,2))

figure(5)
colormap('gray')
image(WH_qp(:,:,3))

figure(6)
colormap('gray')
image(WH_qp(:,:,4))

figure(7)
colormap('gray')
image(WH_qp(:,:,5))


%Nuestro
figure(8)
colormap('gray')
image(WH_Nuestro(:,:,1))

figure(9)
colormap('gray')
image(WH_Nuestro(:,:,2))

figure(10)
colormap('gray')
image(WH_Nuestro(:,:,3))

figure(11)
colormap('gray')
image(WH_Nuestro(:,:,4))

figure(12)
colormap('gray')
image(WH_Nuestro(:,:,5))





