% Equipo: Mariana, Rom�n, Santiago
% 8 de Octubre del 2020
% Proyecto 1 Optimizaci�n
% -------------------------------------------------------------------------
% Obtenci�n de la descomposici�n de im�genes en dos matrices W, H
% v�a el m�todo de descenso en dos pasos a partit del m�todo de puntos
% interiores.
% -------------------------------------------------------------------------
% Configuraci�n inicial
clear; clc; close all; warning('off');

% Obtener imagen en el ambiente local
load('clown')

% obtenci�n de W, H para k = 5, 20, 30, 60, 80
k = [5; 20; 30; 60; 80];

% metodo de optimizaci�n quadprog
% salvar estructuras
WH_qp = zeros(200, 320, 5);
time_qp = zeros(5);

% metodo
for i = 1:numel(k)
    t0 = cputime;
    [W, H] = descenso2pasos_qp(X, k(i));
    WH_qp(:,:,i) = W*H;
    time_cpu = cputime - t0;
    time_qp(i) = time_cpu;
end 

