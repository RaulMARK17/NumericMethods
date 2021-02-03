# Por Raúl A. Góngora Vázquez, estudiante de Ing. Biomédica.
"""
NumericMethods.py(V1) es una libreria que cuenta con varios metodos numericos:

-Interpolación lineal
-Interpolación de LaGrange
-Interpolación de Newton

-Aproximación lineal por mínimos cuadrados

Para localizar raices:
 -Newton-Raphson

-Calculo
 -derivación numérica por diferencia centrada de cuarto orden
 -integración por el metodo del trapecio

-Ecuaciones diferenciales
 -metodo de euler
 -euler mejorado
 -Runge-Kutta de cuarto orden

Se pueden graficar las funciones con los puntos dados, imprimir las funciones resultantes en formato de texto o en codigo.
REQUIERE NUMPY Y MATPLOTLIB
"""

import matplotlib.pyplot as plt
import numpy as np


class laGrangeInterpolation:
    """
    Objeto de interpolación polinomial por LaGrange:\n
    \tLlamar a help() para más información.

    Parametros
    ----------
    \tx = arreglo de puntos en x
    \ty = arreglo de puntos en y
    \tambos de la misma longitud n+1, del cual se obtendra un polinomio de grado n.

    """

    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.steps = len(self.x)

        if (len(self.x) != len(self.y)):
            print('LOS ARREGLOS SON INCORRECTOS.')

    def f(self, x):
        """
        Función resultante
        """
        Px = 0
        for i in range(self.steps):
            u = 1
            v = 1
            for j in range(self.steps):
                if (j != i):
                    u *= (x - self.x[j])
                    v *= (self.x[i] - self.x[j])

            Px += self.y[i]*(u/v)

        return Px

    def getFunct(self):
        """
        Imprime la funcion como texto
        """
        Px = ''
        for i in range(self.steps):
            if ((i > 0) & (i < (self.steps - 1))):
                Px += ' + '
            u = ''
            v = ''
            for j in range(self.steps):
                if (j != i):
                    u += '(x - {})'.format(self.x[j])
                    v += '({} - {})'.format(self.x[i], self.x[j])
            Px += '{}({}/{})'.format(self.y[i], u, v)
        print('P(x) = {}'.format(Px))

    def getShortFunct(self):
        """
        Imprime la funcion como texto, pero corta
        """
        Px = ''
        for i in range(self.steps):
            if ((i > 0) & (i < (self.steps - 1))):
                Px += ' + '
            u = ''
            v = 1
            for j in range(self.steps):
                if (j != i):
                    u += '(x - {})'.format(self.x[j])
                    v *= (self.x[i] - self.x[j])
            Px += '{}({}/{})'.format(self.y[i], u, str(v))
        print('P(x) = {}'.format(Px))

    def getFunctPy(self):
        """
        Imprime la funcion como código
        """
        print('def f(x):')
        Px = ''
        for i in range(self.steps):
            if ((i > 0) & (i < (self.steps))):
                Px += ' + '
            u = ''
            v = 1
            for j in range(self.steps):
                if ((i == 0) & (j == 1)):
                    pass
                else:
                    if ((j > 0) & (j < (self.steps)) & (j != i)):
                        u += '*'
                if (j != i):
                    u += '(x - {})'.format(self.x[j])
                    v *= (self.x[i] - self.x[j])
            Px += '{}*(({})/{})'.format(self.y[i], u, str(v))
        print('    return ({})'.format(Px))

    def graph(self, **kwargs):
        """
        Imprime un grafico con la función y los puntos.

         kwargs(opcionales) 
        -------------------
        \ta = limite a
        \tb = limite b
        \tx = [] (arreglo de nuevos puntos a evaluar, en color morado)
        """
        fix = 10*((np.amax(self.x) - np.amin(self.x))/np.amax(self.x))
        self.a = kwargs.get('a', np.amin(self.x) - fix)
        self.b = kwargs.get('b', np.amax(self.x) + fix)
        self.array = kwargs.get('x', None)

        x = np.linspace(self.a, self.b, 100)
        plt.plot(x, self.f(x))
        for i in range(self.steps):
            plt.plot(self.x[i], self.y[i], 'o', Color='red')
        if self.array:
            for i in range(len(self.array)):
                plt.plot(self.array[i], self.f(
                    self.array[i]), 'o', Color='purple')

        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    def graphPoints(self):
        """
        Imprime un grafico con los puntos.
        """
        for i in range(self.steps):
            plt.plot(self.x[i], self.y[i], 'o', Color='red')

        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    def dataTable(self):
        """
        Imprime una tabla con los puntos incluidos
        """
        z = np.zeros((self.steps, 2))

        for i in range(self.steps):
            z[i][0] = self.x[i]
            z[i][1] = self.y[i]

        print('   x    y\n')
        print(np.matrix(z))


class NewtonInterpolation:
    """
    Objeto de interpolación de Newton:
    \tLlamar a help() para más información.

    Parametros:
    -----------
    \tx = arreglo de puntos en x
    \ty = arreglo de puntos en y
    \tambos deben de ser de la misma longitud

    """

    def __init__(self, x, y):

        self.x = x
        self.y = y
        self.steps = len(self.x)
        self.divsFinitas = []

        if (len(self.x) != len(self.y)):
            print('LOS ARREGLOS SON INCORRECTOS.')

        for i in range(self.steps - 1):
            self.divsFinitas.append([])

        for i in range(self.steps):
            if (i+1 == self.steps):
                break
            else:
                self.divsFinitas[0].append(
                    (self.y[i+1] - self.y[i])/(self.x[i+1] - self.x[i]))

        k = 1
        for i in range(self.steps):
            k += 1
            if(i+1 == self.steps):
                break

            for j in range(self.steps):
                if (j+k == self.steps):
                    break

                else:
                    self.divsFinitas[i+1].append(
                        (self.divsFinitas[i][j+1] - self.divsFinitas[i][j]) / (self.x[j+k] - self.x[j]))

    def f(self, x):
        """
        Función resultante
        """
        y = self.y[0]
        xs = 1
        for i in range(self.steps - 1):
            xs *= (x - self.x[i])
            y += xs*self.divsFinitas[i][0]
        return y

    def getFunct(self):
        """
        Imprime la funcion como texto
        """
        xs = ''
        Str = ''
        for i in range(self.steps - 1):
            xs += '(x - {})'.format(str(self.x[i]))
            Str += '({}){}'.format(str(self.divsFinitas[i][0]), str(xs))
            if (i < (self.steps - 2)):
                Str += ' + '
        print('f(x) = {} + {}'.format(str(self.y[0]), Str))

    def getFunctPy(self):
        """
        Imprime la funcion como código
        """
        print('def f(x):')
        xs = ''
        Str = ''
        for i in range(self.steps - 1):
            if (i > 0):
                xs += '*'
            xs += '(x - {})'.format(str(self.x[i]))
            Str += '({})*{}'.format(str(self.divsFinitas[i][0]), str(xs))
            if (i < (self.steps - 2)):
                Str += ' + '
        print('    return ({} + {})'.format(str(self.y[0]), Str))

    def graph(self, **kwargs):
        """
        Imprime un grafico con la función y los puntos.

         kwargs(opcionales) 
        -------------------
        \ta = limite a
        \tb = limite b
        \tx = [] (arreglo de nuevos puntos a evaluar, en color morado)
        """

        fix = 10*((np.amax(self.x) - np.amin(self.x))/np.amax(self.x))
        self.a = kwargs.get('a', np.amin(self.x) - fix)
        self.b = kwargs.get('b', np.amax(self.x) + fix)
        self.array = kwargs.get('x', None)

        x = np.linspace(self.a, self.b, 100)
        plt.plot(x, self.f(x))
        for i in range(self.steps):
            plt.plot(self.x[i], self.y[i], 'o', Color='red')
        if self.array:
            for i in range(len(self.array)):
                plt.plot(self.array[i], self.f(
                    self.array[i]), 'o', Color='purple')

        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    def graphPoints(self):
        """
        Imprime un grafico con los puntos.
        """
        for i in range(self.steps):
            plt.plot(self.x[i], self.y[i], 'o', Color='red')

        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    def dataTable(self):
        """
        Imprime una tabla con los puntos incluidos
        """
        z = np.zeros((self.steps, 2))

        for i in range(self.steps):
            z[i][0] = self.x[i]
            z[i][1] = self.y[i]

        print('   x    y\n')
        print(np.matrix(z))


class linealInterpolation:
    """
    Objeto de interpolación polinomial Lineal:
    \tLlamar a help() para más información.

    Parametros:
    -----------
    \tx = arreglo de puntos en x
    \ty = arreglo de puntos en y
    \tambos deben de ser de la misma longitud

    """

    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.steps = len(self.x)

        if (len(self.x) != len(self.y)):
            print('LOS ARREGLOS SON INCORRECTOS.')

    def f(self, x):
        """
        La función resultante.

        El valor a ingresar debe estar entre el rango de puntos en x que se entregaron antes o no funcionará.
        """
        self.x1 = 0
        self.x2 = 0
        self.y1 = 0
        self.y2 = 0
        for i in range(self.steps):
            if i+1 == self.steps:
                break

            else:
                if ((x > self.x[i]) and (x < self.x[i+1])):
                    self.x1 = self.x[i]
                    self.x2 = self.x[i+1]
                    self.y1 = self.y[i]
                    self.y2 = self.y[i+1]
                    break

        return ((self.y2 - self.y1)/(self.x2-self.x1))*(x - self.x1) + self.y1

    def getFunct(self):
        """
        Imprime la funcion como texto
        """
        a = '(x - {})'.format(self.x1)
        b = str(self.x2 - self.x1)

        c = str(self.y2 - self.y1)

        print('f(x) = ({}/{}){} + {}'.format(a, b, c, self.y1))

    def getFunctPy(self):
        """
        Imprime la funcion como código
        """
        print('def f(x):')
        a = '(x - {})'.format(self.x1)
        b = str(self.x2 - self.x1)

        c = str(self.y2 - self.y1)

        print('    return ({}/{})*{} + {}'.format(a, b, c, self.y1))

    def graph(self, x=None):
        """
        Imprime un grafico con la función y los puntos.

        Opcional a introducir:
        ---------------------
        \tx = arreglo de nuevos puntos a evaluar, en color morado (entre el rango de los puntos con los que se hizo la interpolación).
        """
        self.a = self.x[0] + .00001
        self.b = self.x[self.steps - 1] - .00001
        self.array = x

        x = np.linspace(self.a, self.b, 100)
        y = []
        for i in x:
            y.append(self.f(i))

        plt.plot(x, y)
        for i in range(self.steps):
            plt.plot(self.x[i], self.y[i], 'o', Color='red')
        if self.array:
            for i in range(len(self.array)):
                plt.plot(self.array[i], self.f(
                    self.array[i]), 'o', Color='purple')

        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    def graphPoints(self):
        """
        Imprime un grafico con los puntos.
        """
        for i in range(self.steps):
            plt.plot(self.x[i], self.y[i], 'o', Color='red')

        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    def dataTable(self):
        """
        Imprime una tabla con los puntos incluidos
        """
        z = np.zeros((self.steps, 2))

        for i in range(self.steps):
            z[i][0] = self.x[i]
            z[i][1] = self.y[i]

        print('   x    y\n')
        print(np.matrix(z))


class msAprox:
    """
    Aproximación lineal por mínimos cuadrados.
    \tLlamar a help() para más información.

    Parametros
    ----------
    \tx = arreglo de puntos en x
    \ty = arreglo de puntos en y
    \tambos deben de ser de la misma longitud
    """

    def __init__(self, x, y):
        sumxy = 0
        sumx = 0
        sumy = 0
        sumx2 = 0
        self.x = x
        self.y = y
        self.steps = len(x)
        if (len(self.x) != len(self.y)):
            print('LOS ARREGLOS SON INCORRECTOS.')

        n = len(x)
        for i in range(0, len(x)):
            sumxy += x[i]*y[i]
            sumx += x[i]
            sumy += y[i]
            sumx2 += x[i]**2
        self.a = (n*sumxy - sumx*sumy)/(n*sumx2 - sumx**2)
        self.b = (sumy - self.a*sumx)/n

    def f(self, x):
        """
        Función de la recta resultante
        """
        return self.a*x + self.b

    def graphPoints(self):
        """
        Imprime un grafico con los puntos.
        """

        for i in range(self.steps):
            plt.plot(self.x[i], self.y[i], 'o', Color='red')

        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        plt.show()

    def graph(self, **kwargs):
        """
        Imprime una tabla con los puntos incluidos

        kwargs(opcionales) 
        -------------------
        \ta = limite a
        \tb = limite b
        \tx = [] (arreglo de nuevos puntos a evaluar, en color morado)
        """

        fix = 2*((np.amax(self.x) - np.amin(self.x))/np.amax(self.x))
        a = kwargs.get('a', np.amin(self.x) - fix)
        b = kwargs.get('b', np.amax(self.x) + fix)
        self.array = kwargs.get('x', None)

        x = np.linspace(a, b, 100)
        plt.plot(x, self.f(x))
        for i in range(self.steps):
            plt.plot(self.x[i], self.y[i], 'o', Color='red')
        if self.array:
            for i in range(len(self.array)):
                plt.plot(self.array[i], self.f(
                    self.array[i]), 'o', Color='purple')

        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis("equal")
        plt.show()

    def getFunct(self):
        """
        Imprime la funcion como texto
        """
        if self.b > 0:
            print('f(x) = {}x + {}'.format(self.a, self.b))
        else:
            print('f(x) = {}x - {}'.format(self.a, abs(self.b)))

    def getFunctPy(self):
        """
        Imprime la funcion como código
        """
        print('def f(x):')

        if self.b > 0:
            print('    return {}*x + {}'.format(self.a, self.b))
        else:
            print('    return {}*x - {}'.format(self.a, abs(self.b)))

    def dataTable(self):
        """
        Imprime una tabla con los puntos incluidos
        """
        z = np.zeros((self.steps, 2))

        for i in range(self.steps):
            z[i][0] = self.x[i]
            z[i][1] = self.y[i]

        print('   x    y\n')
        print(np.matrix(z))


class euler:
    def __init__(self, ED, x0, y0, **kwargs):
        """
        Metodo de euler, entrega un codigo para copiar y ejecutar.
        Llamar a help() para más información.

        Parametros
        ----------
        \tED = funcion en string en formato python
        \tx0 = punto inicial x
        \ty0 = punto inicial y


        kwargs(opcionales) 
        -------------------
            \th = tamaño de paso
            \tn = número de pasos(por defecto 10000)
            \txn = x final
        """
        self.h = kwargs.get('h', None)
        self.n = kwargs.get('n', 10000)
        xn = kwargs.get('xn', None)

        if ((self.h != None) and (xn != None)):
            self.n = round((xn - x0)/self.h)

        if ((self.h == None) and (xn != None)):
            self.h = (xn - x0)/self.n

        self.x0 = x0
        self.y0 = y0
        self.ED = ED

        print(
            """
import matplotlib.pyplot as plt

def ED(x, y):
    return {}

n = {}
x0 = {}
y0 = {}
h = {}
y_ = [y0]
x = [x0]

for i in range(n):
    yi = y_[i] + ED(x[i], y_[i])*h
    y_.append(yi)
    x.append(x[i] + h)
print('x = {}, y = {}'.format(x[-1], y_[-1]))

plt.plot(x, y_)
plt.show()    
        """.format(self.ED, self.n, self.x0, self.y0, self.h, "{}", "{}"))

    def compare(self, f):
        """
        Si tienes la Ecuación diferencial resuelta:

        Parametros
        --------------
        \tf = funcion en string, como codigo python
        """
        print(
            """
def f(x):
    return {}

plt.plot(x, y_)
plt.plot(x, f(x), c='red')
plt.show()

print('Error = {}%'.format((f(x[-1]) - y_[-1])/f(x[-1])))
            """.format(f, '{}'))


class eulerImproved:
    def __init__(self, ED, x0, y0, **kwargs):
        """
        Metodo de euler mejorado, entrega un codigo para copiar y ejecutar.
        Llamar a help() para más información.

        Parametros
        ----------
        \tED = funcion en string en formato python
        \tx0 = punto inicial x
        \ty0 = punto inicial y


        kwargs(opcionales) 
        -------------------
            \th = tamaño de paso
            \tn = número de pasos(por defecto 10000)
            \txn = x final
        """
        self.h = kwargs.get('h', None)
        self.n = kwargs.get('n', 10000)
        xn = kwargs.get('xn', None)

        if ((self.h != None) and (xn != None)):
            self.n = round((xn - x0)/self.h)

        if ((self.h == None) and (xn != None)):
            self.h = (xn - x0)/self.n

        self.x0 = x0
        self.y0 = y0
        self.ED = ED

        print(
            """       
import matplotlib.pyplot as plt

def ED(x, y):
    return {}

n = {}
x0 = {}
y0 = {}
h = {}
y_ = [y0]
x = [x0]

for i in range(0, n):
    x.append(x[i] + h)
    yi_ = y_[i] + h*ED(x[i], y_[i])
    yi = y_[i] + h*(ED(x[i], y_[i]) + ED(x[i+1], yi_))/2
    y_.append(yi)
    
print('x = {}, y = {}'.format(x[-1], y_[-1])) 

plt.plot(x, y_)
plt.show()    
        """.format(self.ED, self.n, self.x0, self.y0, self.h, "{}", "{}"))

    def compare(self, f):
        """
        Si tienes la Ecuación diferencial resuelta:

        Parametros
        --------------
        \tf = funcion en string, como codigo python
        """
        print(
            """
def f(x):
    return {}

plt.plot(x, y_)
plt.plot(x, f(x), c='red')
plt.show()

print('Error = {}%'.format((f(x[-1]) - y_[-1])/f(x[-1])))
            """.format(f, '{}'))


class rungeKutta_4:
    def __init__(self, ED, x0, y0, **kwargs):
        """
        Runge-Kutta de cuarto orden, entrega un codigo para copiar y ejecutar.
        Llamar a help() para más información.

        Parametros
        ----------
        \tED = funcion en string en formato python
        \tx0 = punto inicial x
        \ty0 = punto inicial y


        kwargs(opcionales) 
        -------------------
            \th = tamaño de paso
            \tn = número de pasos(por defecto 10000)
            \txn = x final
        """
        self.h = kwargs.get('h', None)
        self.n = kwargs.get('n', 10000)
        xn = kwargs.get('xn', None)

        if ((self.h != None) and (xn != None)):
            self.n = round((xn - x0)/self.h)

        if ((self.h == None) and (xn != None)):
            self.h = (xn - x0)/self.n

        self.x0 = x0
        self.y0 = y0
        self.ED = ED

        print(
            """
import matplotlib.pyplot as plt

def ED(x, y):
    return {}
    
n = {}
x = [{}]
y_ = [{}]
h = {}

for i in range(n):
    k1 = h*ED(x[i], y_[i])
    k2 = h*ED(x[i] + 0.5*h, y_[i] + 0.5*k1)
    k3 = h*ED(x[i] + 0.5*h, y_[i] + 0.5*k2)
    k4 = h*ED(x[i] + h, y_[i] + k3)

    y_.append(y_[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4))
    x.append(x[i] + h)


plt.plot(x, y_)
plt.show()
                """.format(self.ED, self.n, self.x0, self.y0, self.h))

    def compare(self, f):
        """
        Si tienes la Ecuación diferencial resuelta:

        Parametros
        --------------
        \tf = funcion en string, como codigo python
        """
        print(
            """
def f(x):
    return {}

plt.plot(x, y_)
plt.plot(x, f(x), c='red')
plt.show()

print('Error = {}%'.format((f(x[-1]) - y_[-1])/f(x[-1])))
            """.format(f, '{}'))


class calculus:
    """
    Métodos de integracion y diferenciacion numericos:\n
    \tLlamar a help() para más información.
    """

    def integrate(self, f, a, b, **kwargs):
        """
        Entrega un codigo para copiar y ejecutar.

        Parametros
        ----------
        \tf = funcion en string en formato python
        \ta = incio de intervalo
        \tb = final de intervalo

        kwargs(opcionales) 
        -------------------
            \th = tamaño de rectangulo
            \tn = número de rectangulos(por defecto 10000)
        """
        h = kwargs.get('h', None)
        n = kwargs.get('n', 10000)

        if ((h != None) and (n == 10000)):
            n = (b - a)/h

        print(
            """
a = {}
b = {}
n = {}

def f(x):
    return {}

def integrate(a, b, n):
    _ = (b-a)/n
    increase = a
    integral = 0
    for i in range(n):
        integral += (f(increase) + f(increase + _))*((increase + _) - increase)/2
        increase += _
    return integral

print(integrate(a, b, n))
        """.format(a, b, n, f))

    def derive(self, f, x, h):
        """
        Entrega un codigo para copiar y ejecutar.

        Parametros
        ----------
        \tf = funcion en string en formato python
        \tx = punto a evaluar
        \th = aproximación a 0
         """
        print(
        """
x = {}
h = {}
def f(x):
    return {}

def derive(x, h):
    return (-f(x + 2*h) + 8*f(x+h) - 8*f(x-h) + f(x - 2*h))/(12*h)

print(derive(x, h))
            """.format(x, h, f))


class newtonRaphson:
    def __init__(self, x0, f, **kwargs):
        """
        Entrega un codigo para copiar y ejecutar.

        Parametros
        ----------
        \tf = funcion en string en formato python
        \tx0 = punto de inicio

        kwargs(opcionales) 
        -------------------
            \tderivative = función string derivada en formato de codigo python(si no se proporciona, se aproximara numericamente)
            \tlimit = limite de iteraciones( por defecto 10000)
            \tp = precisión del método(por defecto 0.00001)
        """
        derivative = kwargs.get('derivative', None)
        limit = kwargs.get('limit', 10000)
        p = kwargs.get('p', 0.00001)

        if derivative == None:
            print(
                """
import matplotlib.pyplot as plt
import numpy as np

x = {}
limit = {}
p = {}

def f(x):
    return {}

def derivative(x):
    h = 0.000001
    return (-f(x + 2*h) + 8*f(x+h) - 8*f(x-h) + f(x - 2*h))/(12*h)

iteration = 0
while((f(x) > p) or (f(x) < -p) or (x < 0)):
    xi = (x - f(x)/derivative(x))
    x = xi
    iteration += 1
    if iteration >= limit:
        break

print('x = {}'.format(x))
print('f(x) = {}'.format(f(x)))
print('iteración: {}'.format(iteration))

X = np.linspace(-1+x, x+1, 100)

plt.plot(X, f(X))
plt.plot(x, f(x), 'o', c='red')
plt.show()
            """.format(x0, limit, p, f, '{}', '{}', '{}'))

        else:
            print(
                """
import matplotlib.pyplot as plt
import numpy as np

x = {}
limit = {}
p = {}

def f(x):
    return {}

def derivative(x):
    return {}

iteration = 0
while((f(x) > p) or (f(x) < -p) or (x < 0)):
    xi = (x - f(x)/derivative(x))
    x = xi
    iteration += 1
    if iteration >= limit:
        break

print('x = {}'.format(x))
print('f(x) = {}'.format(f(x)))
print('iteración: {}'.format(iteration))

X = np.linspace(-1+x, x+1, 100)

plt.plot(X, f(X))
plt.plot(x, f(x), 'o', c='red')
plt.show()
            """.format(x0, limit, p, f, derivative, '{}', '{}', '{}'))


def help():
    print(
    """
    Las clases de interpolación y msAprox tienen las siguientes funciónes:\n
    \tgetFunct(): obtener la función en texto.
    \tgetFunctPy(): obtener la función en código python.
    \tgraph(): Gráfica con los puntos y la función
    \tgraphPoints(): Gráfica con los puntos dados.
    \tf(): Para llamar a la función resultante y evaluar culaquier valor de x.
    \tdataTable(): Para ver en formatos de tabla los puntos introducidos(msAprox no la tiene)

    La clase de calculo entrega codigos para copiar, pegar y ejecutar; tiene las siguientes funciones:\n
    \tintegrate(): requiere de la función, un punto y los limites para encontrar la integral definida.
    \tderive(): requiere de la función, un punto y la h para encontrar la derivada en el punto dado.

    La clase newtonRaphson entrega un código para copiar, pegar y ejecutar; requiere:\n
    \tf: la función.
    \tx0: el punto inicial.
    \topcionalmente, se puede ajustar la precisión e introducir la función derivada.

    Las clases de ecuaciones diferenciales entregan códigos para copiar, pegar y ejecutar; requiere:\n
    \tED: la ecuacion diferencial en terminos de x e y.
    \tx0: la coordenada x inicial.
    \ty0: la coordenada y inicial.
    \topcionalmente, se puede ajustar el tamaño de paso, numero de pasos y la x final
    \t también tiene una funcion para comparar la grafica aproximada con la real si se cuenta con la funcion resultado de la ED.

    """)
