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
% Optimizaci�n Num�rica
% Proyecto 1
% 30 de septiembre del 2020
% ITAM
%--------------------------------------------------------------------------
[r,p] = size(X);    % obtenemos las dimensiones de X
W = ones(r,k);      % inicializamos a W
H = ones(k,p);       % inicializamos a H
maxiter = 3;
for kiter = 1:maxiter     % corremos el m�todo k veces
   % El probelma de min ||X - Wk*H||_F^2 equivale a p problemas cuadr�ticos
   % de la forma: min H*j'(Wk'Wk)H*j - X*j'WkH*j s.a. H*j >= 0 para j=1,...,p
   for j = 1:p
       columnasX = X(:,j);
       Q1 = (W'*W); %dimension k*k
       c1 = (-columnasX'*W)'; % dimension k*1
       A1 = eye(k);
       b1 = zeros(k,1);
       % Usando el metodo de puntos interiores
       [x1,y1,mu1] = punintpc(Q1, A1, c1, b1);
       % sol. es H = [x1|x2|...|xp]
       H(:,j)= x1;
   end
   % El probelma de min ||X - W*Hk||_F^2 equivale a r problemas cuadr�ticos
   % de la forma: min Wi*(Hk*Hk')Wi*' - Xi*Hk'Wi*' s.a. Wi* >= 0 para i=1,...,r
   for i = 1:r
       renglonesX = X(i, :); %dimension 1*p
       Q2 = (H*H'); %dimension k*k
       c2 = (-renglonesX*H')'; % dimension k*1
       A2 = eye(k);
       b2 = zeros(k,1);
       % Usando el metodo de puntos interiores
       [x2,y1,mu2] = punintpc(Q2, A2, c2, b2);
       % sol. es W' = [x1|x2|...|xr]
       W(i, :) = x2';
   end
end    



