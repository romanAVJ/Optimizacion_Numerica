function [W, H] = descenso2pasos(X, k)
% Metodo de descenso en dos pasos 
% Min   ||X - W*H ||_F^2
% s.a.   W >= 0
%        H >= 0
%
% Llamado: [W, H] = descenso2pasos(X, k);
% Input
% X.- matriz rxp con X(i,j) >= 0
% k.- numero de iteraciones del metodo
%
%Out
% W.- matriz rxk con W(i,j) >= 0
% H.- matriz kxp con H(i,j) >= 0
% con W*H aproxima la matriz X
%
%--------------------------------------------------------------------------
% Optimización Numérica
% Proyecto 1
% 30 de septiembre del 2020
% ITAM
%--------------------------------------------------------------------------
[r,p] = size(X);
W = ones(r,k);
H = one(k,p);
for kiter = 1:k
   for j = 1:p
       Q1 = 2*W'.*W;
       c1 = -X(1:r,j)'.*W;
       A1 = eye(k);
       b1 = zeros(k,1);
       % Usando el metodo de puntos interiores
       [x1,~,~] = punintpc(Q1, A1, c1, b1);
       H(1:k,j)= x1;
   end
   
   for i = 1:r
       Q2 = 2*H.*H';
       c2 = -X(i,1:p).*H';
       A2 = eye(k);
       b2 = zeros(k,1);
       % Usaremos el metodo de puntos interiores
       [x2,~,~] = punintpc(Q2, A2, c2, b2);
       W(i,1:k) = x2';
   end
     
end    



