# NumericMethods
Esta librería es una recopilación de códigos que realicé en mi tiempo estudiado métodos numéricos y aplicándolos en algunas otras materias. 
Todos los códigos son automáticos, practicamente cualquier problema de aproximación, derivación, integración es posible resolver mientras esté en los límites y naturaleza del método. 

Pueden graficarse las funciones, los puntos, evaluar nuevos puntos y ajustar los límites de las gráficas.
Obtener las funciones resultantes(de los métodos de interpolación) en formato de código python o texto.

Para ecuaciones diferenciales es posible obtener una gráfica comparativa de la función aproximada y la funcíon resultado de la ecuación diferencial, con un error estimado; siguiendo la formula: (Y - y)/Y, donde Y es el valor real e y es el valor estimado. 

Con Newton-Raphson es posible visualizar que efectivamente el punto ha llegado a donde se esperaba por medio de una gráfica.

Se pueden ajustar límites, precisión, número de pasos, tamaño de pasos, entre otras opciones. 

Todo está comentado y es bastante sencillo de usar, con solo ingresar unos cuantos datos ya se obtienen resultados.

Los métodos presentes son:
-Interpolación lineal
-Interpolación de LaGrange
-Interpolación de Newton

-Aproximación lineal por mínimos cuadrados

-Newton-Raphson

-Calculo
 -derivación numérica por diferencia centrada de cuarto orden
 -integración por el metodo del trapecio

-Ecuaciones diferenciales:
-metodo de euler
-euler mejorado
-Runge-Kutta de cuarto orden
 
 Cada método es una clase. 
 Por la naturaleza de los métodos y las limitaciones del lenguaje(o mi conocimiento de él) algunas clases imprimen un código ejecutable, listo para copiar y ejecutar donde se desee. 
 REQUIERE NUMPY Y MATPLOTLIB.
