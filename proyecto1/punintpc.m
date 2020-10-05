function [x, y, mu] = punintpc(Q,A,c,b)
% Metodo de punto interior para el problema cuadratico
% Min   (0.5)* x' * Q * x + c'* x
% s.a.   A * x >= b
%
%  Llamado: [x, y, mu] = punintpc(Q,A,c,b);
% In
% Q.- matriz nxn simetrica y positiva definida
% A.- matriz mxn con m <= n y rango(A) = m.
% b.- vector columna en R^m .
% c.- vector columna en R^n .
%
%Out
% x.- vector en R^n con la aproximacion del minimo local.
% mu.- vector en R^m con la aproximacion al multiplicador
%          de Lagrange asociado a la restriccion de desigualdad.
% y.- vector en R^m con la variable de holgura en la restriccion
%     de desigualdad.
%---------------------------------------------------------------------------
% Optimizacion Numerica
% Ultima actualizacion 1.oct.20
% Equipo: Santiago Muriel
%         Mariana G Martinez
%         Roman Velez
% ITAM
%--------------------------------------------------------------------------------
% Parametros iniciales
tol = 1e-06;       % Tolerancia a las condiciones necesarias de 1er orden
maxiter = 250;     % maximo numero de iteraciones permitidas
iter = 0;          % contador de las iteraciones
%-----------------------------------------------------------
n = length(c);     % dimension de la variable principal
m = length(b);     % numero de restricciones de desigualdad
%----------------------------------------------------------
% variables iniciales
% Teniamos varias opciones, hemos decidido iniciarlas todas
% en los vectores de unos por facilidad y asi ahorrarnos
% el problema de minimizacion
x = ones(n,1); % minimo aproximado inicial
mu = ones(m,1); % multiplicador de Lagrange inicial
y = ones(m,1); % variable de holgura inicial
gamma = (0.5)*(mu'*y)/m; % actualizacion de la perturbacion de las KKT
%-----------------------------------------------
% vectores para graficacion
cnpo=[]; comp =[];
% Norma de las condiciones necesarias de primer orden
% Como no tenemos restricciones de igualdad, esta H se ve 
% un poco distinta (no tiene un bloque de renglon)
H = [Q*x - A'*mu + c; A*x - y - b; mu.*y];
norma = norm(H);
disp('Iter      CNPO             gamma ')
disp('-----------------------------------------')
while(norma > tol & iter < maxiter) % Parte iterativa del método de Newton
  % Resuelve el sistema lineal de Newton para la trayectoria central
  YminU = (1./y).*mu;
  rx = Q*x - A'*mu + c;
  ry = A*x - y - b;
  rmu = y.*mu - gamma;
  %---------------------------------------------------------- 
    % Resolvemos el sistema lineal
  K = Q + A'*diag(YminU)*A;
  ld = -(rx + A'*diag(YminU)*ry + A'*((1./y).*rmu));
  Dx = K \ ld;
  
  %---------------------------------------------------------- 
    % Sustituimos para los demas pasos
    
    Dy = A*Dx + ry;
    Dmu = -(1./y).*(mu.*Dy + rmu);
  
  
   %---------------------------------------------------------- 
    % Acorta el paso
    bt = []; gm = [];
    for k =1:m
        if (Dmu(k) < 0)
            bt = [bt; -(mu(k)/Dmu(k))];
        else
            bt = [bt; 1];
        end
        if(Dy(k) < 0)
            gm = [gm; -(y(k)/Dy(k))];
        else
            gm = [gm; 1];
        end
    end
    
    alfa = min([bt ; gm]); %alfa que usamos para recortar el paso de las variables
    alfa =(0.9995)*min([1 alfa]);  % garantiza que las mu y z siempre sean positivas
    %-----------------------------------------------------------
     % Nuevo punto
       x      = x + alfa*Dx;
       mu     = mu + alfa*Dmu;
       y      = y + alfa*Dy;
     %-------------------------------------------------------  
     % Nueva tau
        gamma = (0.5)*(mu'*y)/m; % promedio de la complementaridad
     %-------------------------------------------------------  
       %Condiciones necesarias de primer orden
       H = [Q*x - A'*mu + c; A*x - y - b; mu.*y];
       norma = norm(H);
       iter = iter + 1;
       cnpo =[cnpo norma];
       comp = [comp 2*gamma];
       disp(sprintf('%3.0f  %2.8f  %2.8f',iter,norma,2*gamma))
end

   semilogy([1:iter],cnpo,'r',[1:iter],comp,'b')
   title('Convergencia de puntos interiores')
   legend('CNPO', 'Complementaridad')
        
