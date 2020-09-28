
# Optimización Numérica
Proyectos de la materia de optimización numérica impartida por Zeferino Parada, Agosto-Diciemrbe 2020. 

# Equipo
El equipo está conformado por 
- Santiago Muriel
- Mariana Martínez
- Román Vélez

## Actualización 28.sep.20
Se agregó el archivo qpintpoint.m el cual es el Metodo de punto interior para el problema cuadrático:
  Min   (0.5)* x' * Q * x + c'* x
  s.a.   A * x = b
         F*x >= d
         
Debemos de modificar este método, pues lo que queremos resolver es:
  Min   (0.5)* x' * Q * x + c'* x
  s.a.   A * x >= b
  
Por esta modificación, nuestra matriz H (matriz definida por bloques) tiene otra forma:
H = [Q 0 -A' ; A -Im 0 ; 0 U Y]
Por ende la parte iterativa del paso de Newton también tiene otra forma (los cambios están
en la página 2 del pdf de la descipción del proyecto).
Además tenemos que implementar este código 2 veces para el método de descenso en dos pasos.


## Qué falta hacer

 - Modificar el archivo qpintpoint.m, es decir cambiar las restricciones, la matriz H, el paso de Newton
  y la sigma
 - Hacer la función descenso2pasos.m en la cual se implementa el metodo de puntos interiores
 - Hacer un script con la implementación del método de descenso por 2 pasos aplicado a la foto
  del payaso
