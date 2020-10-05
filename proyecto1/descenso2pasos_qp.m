function [W, H] = descenso2pasos_qp(X, k)
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
[r,p] = size(X);    % obtenemos las dimensiones de X
W = ones(r,k);      % inicializamos a W
H = ones(k,p);       % inicializamos a H
for kiter = 1:k     % corremos el método k veces
   % El probelma de min ||X - Wk*H||_F^2 equivale a p problemas cuadráticos
   % de la forma: min H*j'(Wk'Wk)H*j - X*j'WkH*j s.a. H*j >= 0 para j=1,...,p
   for j = 1:p
       Q1 = 2*(W'*W);
       c1 = -X(1:r,j)'*W;
       A1 = eye(k);
       b1 = zeros(k,1);
       % Usando el metodo de puntos interiores
       x1 = quadprog(Q1, c1, -A1, -b1);
       % sol. es H = [x1|x2|...|xp]
       H(1:k,j)= x1;
   end
   % El probelma de min ||X - W*Hk||_F^2 equivale a r problemas cuadráticos
   % de la forma: min Wi*(Hk*Hk')Wi*' - Xi*Hk'Wi*' s.a. Wi* >= 0 para i=1,...,r
   for i = 1:r
       Q2 = 2*(H*H');
       c2 = -X(i,1:p)*H';
       A2 = eye(k);
       b2 = zeros(k,1);
       % Usando el metodo de puntos interiores
       x2 = quadprog(Q2, c2, -A2, -b2);
       % sol. es W' = [x1|x2|...|xr]
       W(i,1:k) = x2';
   end
end    
