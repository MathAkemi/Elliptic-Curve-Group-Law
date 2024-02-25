import matplotlib.pyplot as plt
import numpy as np
from sympy import var, plot_implicit, lambdify


################################
### ELLIPTIC CURVE GROUP LAW ###
################################

#######################################################################################################
# This program is a full implementation of the group law on elliptic curves with plotting. The goal   #
# of writing it is to deepen my understanding of elliptic curves and to improve my ability to work    #
# with them to prepare for further work with elliptic curve cryptography. I intend forthis code to be #
# useful to me as I write implementations for more advanced algorithms.				                  #
#######################################################################################################



#############################################################################################
# The Point class holds a point on the curve.                                               #
#											                                                #
# Point.x: the x-coordinate of the point						                            #
# Point.y: the y-coordinate of the point						                            #
# Point.curve: the EllipticCurve object corresponding to the curve on which the point lies. #
#											                                                #
# Point addition, negation, subtraction, and mult. as a groupf action of the integers are   #
# all fully supported with the respective magic methods.                                    #
#############################################################################################
class Point():
    
    def __init__(self, x, y, E=None):
        if isinstance(E, EllipticCurve):
            self.curve = E
        self.x = x
        self.y = y
    
    def __neg__(self):
        if self.x == float("inf") or self.y == float("inf"):
            return Point(self.x, self.y, self.curve)
        return Point(self.x, -self.y, self.curve)

    def __add__(self, Q):
        return self.curve.addPoints(self, Q)

    def __sub__(self, Q):
        return self.curve.addPoints(self, -Q)

    def __mul__(self, n):
        if isinstance(n, int):
            return self.curve.mulPoint(self, n)

    __lmul__ = __mul__
    __rmul__ = __mul__

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __str__(self):
        return str("(" + str(self.x) + ", " + str(self.y) + ")") 

    __repr__ = __str__



##################################################################################################
# The EllipticCurve class represents the elliptic curve to be studied.                           #
#												                                                 #
# EllipticCurve.a: the coefficient of the x^3 term in the Weierstrass equation                   #
# EllipticCurve.b: the coefficient of the x term in the Weierstrass equation                     #
# EllipticCurve.dis: the discriminant of the curve, must be != 0                                 #
# EllipticCurve.eqn: the equation to output: 'w' for Weierstrass, 'l' for Legendre               #
# EllipticCurve.l: the value of lambda in the Legendre equation                                  #
#                                                                                                #
# isOnCurve (P): returns 1 if P is a point on the curve, 0 otherwise                             #
# addPoints(P, Q): returns P + Q according to the group law on the curve                         #
# timesTwo(P): returns 2*P = P + P								                                 #
# mulPoint(n, P): adds P to itself n-1 times and returns the result using an efficient algorithm #
# plotCurve(): plots the elliptic curve on the real plane                                        #
# addAndPlot(P, Q): computes P + Q and plots the process of finding the value			         #
#												                                                 #
# Comparison with "==" is supported. Calling the object returns the Weierstrass equation.        #
##################################################################################################
class EllipticCurve():

    # Initialize an EllipticCurve object.
    def __init__(self, a, b, eq='w'):
        self.a = a  # Weierstrass constants; eqn type only specifies the type to print
        self.b = b
        if eq == 'w' or eq == 'l' or eq == 'c' or eq == 'q':
            self.eqn = eq
        else:
            print("Equation type invalid. Please choose either 'w' for Weierstrass, 'l' for Legendre, 'c' for cubic, and 'q' for quartic equations. Set to Weierstrass by default.")
            self.eqn = 'w'
        self.l = 0

        self.dis = -16*(4*a*a*a+27*b*b) 
        if self.dis == 0:
            raise Exception("This is not an elliptic curve.")

    # Reutrn 1 if P is on the curve and 0 otherwise.
    def isOnCurve(self, P):
        toCheck = (P.y)**2 - (P.x)**3 - self.a*(P.x) - self.b
        if P.x == float("inf") or P.y == float("inf"):
            return 1 
        if round(toCheck, 5) == 0:
            return 1
        else:
            return 0

    # Add the points P and Q according to the group law on the curve E.
    def addPoints(self, P, Q):
        if not isinstance(P, Point) or not isinstance(Q, Point):
            return "Please enter a Point object."
        if not self.isOnCurve(P) or not self.isOnCurve(Q):
            return "Please enter a Point on the curve."

        if not P.x == float("inf") and not P.y == float("inf") and not Q.x == float("inf") and not Q.y == float("inf") and not (P.x == Q.x):
            m = (Q.y - P.y) / (Q.x - P.x)
            x = m*m - Q.x - P.x
            y = -(m*(x - P.x) + P.y)

        if P.x == Q.x:
            if P.y != Q.y:
                return Point(float("inf"), float("inf"), curve)
            else:
                return curve.timesTwo(P)

        if P.x == float("inf") and P.y == float("inf"):
            x = Q.x
            y = Q.y

        if Q.x == float("inf") and Q.y == float("inf"):
            x = P.x
            y = P.y

        R = Point(x, y, curve)
        return R

    # Double a point, i.e. add it to itself.
    def timesTwo(self, P):
        m = (3 * P.x * P.x + P.curve.a) / (2 * P.y)
        
        x = m*m - 2*P.x
        y = -(m*(x - P.x) + P.y)
        Q = Point(x, y, curve)

        return Q

    # Multiply a point by an integer n, i.e. add it to itself n-1 times.
    def mulPoint(self, P, n):
        if P.x == float("inf") or P.y == float("inf"):
            return P
        
        a = n
        B = Point(float("inf"), float("inf"), curve)
        C = P

        while a != 0:
            if a % 2 == 0:
                a = a/2
                B = B
                C = self.timesTwo(C)
            elif a % 2 == 1:
                a = a-1
                B = B + C
                C = C

        return B
    
    # Plot the curve.
    def plotCurve(self):
            var('x y')
            toPlot = x**3 + self.a*x + self.b - y**2 
            plot_implicit(toPlot)

    # Add P and Q and plot them on the curve.
    def addAndPlot(self, P, Q):
        R = P + Q

        if not P.x == float("inf") and not P.y == float("inf") and not Q.x == float("inf") and not Q.y == float("inf") and not (P.x == Q.x):
            m = (Q.y - P.y) / (Q.x - P.x)
            c = P.y - m*P.x

        if P.x == Q.x:
            if P.y != Q.y:
                m = float("inf")
            else:
                m = (3 * P.x * P.x + P.curve.a) / (2 * P.y)
                c = P.y - m*P.x

        case = -1   # This is used to solve problems with addition w/ inf.
        if P.x == float("inf") and P.y == float("inf"):
            case = 0
            m = float("inf")

        if Q.x == float("inf") and Q.y == float("inf"):
            case = 1
            m = float("inf")

        var('x y')
        if m != float("inf"):
            x_values = np.linspace(min(P.x, R.x, Q.x) - 1, max(P.x, R.x, Q.x) + 1, 10000)
            y_values = m*x_values + c
        elif m == float("inf"):
            x_values = np.linspace(R.x - 3, R.x + 3, 10000)

        f_curve = lambda x : x**3 + curve.a*x + curve.b
        curve_values_pos = np.sqrt(f_curve(x_values))
        curve_values_neg = -np.sqrt(f_curve(x_values))    
        plt.plot(x_values, curve_values_pos, color='orange')
        plt.plot(x_values, curve_values_neg, color='orange')

        if m != float("inf"):
            plt.plot(x_values, y_values)
        else:
            plt.axvline(x=R.x)
        plt.axvline(x=R.x, linestyle='--', color='gray')
        plt.scatter([P.x, Q.x, R.x], [P.y, Q.y, R.y], color=['red', 'green', 'blue'])
        
        plt.axis('equal')
        plt.show()

    def printEquation(self, c):
        if not isinstance(c, char) or not (c == 'w' or c == 'l'):
            return "Please pass either 'w' or 'l' for the Weierstrass or Legendre equations."
        if c == 'w':
            return __str__(self)
        if c == 'l':
            l = 0;       # WIP: will print the Legendre form of the curve.
            return l;

    def __eq__(self, other):
        return self.a == other.a and self.b == other.b

    def __str__(self):
        if self.a > 0 and self.b > 0:
            return str("y^2 = x^3 + " + str(self.a) + "x + " + str(self.b))
        elif self.a > 0 and self.b < 0:
            return str("y^2 = x^3 + " + str(self.a) + "x - " + str(-self.b))
        elif self.a < 0 and self.b > 0:
            return str("y^2 = x^3 - " + str(-self.a) + "x + " + str(self.b))
        elif self.a < 0 and self.b < 0:
            return str("y^2 = x^3 - " + str(-self.a) + "x - " + str(-self.b))

    __repr__ = __str__


curve = EllipticCurve(-1, 1)
P = Point(-1, 1, curve)
Q = Point(0, 1, curve)
inf = Point(float("inf"), float("inf"), curve)
