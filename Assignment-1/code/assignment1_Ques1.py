from sympy import *
from numpy import *
import math
import matplotlib.pyplot as plt

def plotFuntionClosed(eq,l,u):
    x=Symbol('x')
    f=lambdify(x,eq,"numpy")
    x_vals = linspace(l,u, 100)
    y_vals = f(x_vals)
    plt.plot(x_vals, y_vals)
    plt.title('Plot of F(x) vs x')
    plt.legend(['y='+ eq])
    plt.grid()
    plt.show()
    return 0

def plotFuntionOpen(eq):
    x=Symbol('x')
    f=lambdify(x,eq,"numpy")
    x_vals = linspace(float(input('Enter the lower value of range of graph on x axis: ')),float(input('Enter the upper value of range of graph on x axis: ')), 100)
    y_vals = f(x_vals)
    plt.plot(x_vals, y_vals)
    plt.title('Plot of F(x) vs x')
    plt.legend(['y='+ eq])
    plt.grid()
    plt.show()
    return 0

def bisection(equation):
    x=Symbol('x')
    eq=sympify(equation)
    lower=float(input('Please enter lower starting value (eg. 0.1)\n'))
    upper=float(input('Please enter upper starting value (eg. 1)\n'))
    lo=lower #declaring lo and up for use in 'plotFuntion' funtion.
    up=upper
    f_lower=eq.subs(x,lower).evalf()
    f_upper=eq.subs(x,upper).evalf()
    if f_lower*f_upper>0:
        print("Error: Please select correct interval")
        quit()
    maxiters = int(input('Enter the Maximum Iterations:'))
    maxerror = float(input('Enter the Maximum Error (%):'))
    errors=[math.nan]
    x0=lower #Previous value of xr.
    for i in range(0,maxiters):
        xr=(lower+upper)/2 #updating the current value of xr
        f_xr=eq.subs(x,xr).evalf()
        a=f_lower*f_xr
        if i:
            err=abs(((xr-x0)/xr)*100) #Relative percentage error
            errors.append(err)
        else:
            errors.append(math.nan)
        if a==0:
            break
        elif a>0:
            lower=xr
            f_lower=f_xr
        else:
            upper=xr
            f_upper=f_xr
        if i:
            if (err<maxerror or f_xr==0):
                break
        x0=xr
    print(f"The root of the given funtion is {xr} ")
    plt.plot(errors,marker='o')
    plt.xlabel('Iteration number')
    plt.ylabel('Relative Error(%)')
    plt.title('Plot of relative errors vs the iteration number.')
    plt.grid()
    plt.show()
    plotFuntionClosed(equation,lo,up)
    return 0

def false_position(equation):
    x=Symbol('x')
    eq=sympify(equation)
    lower=float(input('Please enter lower starting value (eg. 0.1)\n'))
    upper=float(input('Please enter upper starting value (eg. 1)\n'))
    lo=lower #declaring lo and up for use in 'plotFuntion' funtion.
    up=upper
    f_lower=eq.subs(x,lower).evalf()
    f_upper=eq.subs(x,upper).evalf()
    if f_lower*f_upper>0:
        print("Error: Please select correct interval")
        quit()
    maxiters = int(input('Enter the Maximum Iterations:'))
    maxerror = float(input('Enter the Maximum Error (%):'))
    errors=[math.nan]
    x0=lower
    for i in range(0,maxiters):
        xr=(f_lower*upper-f_upper*lower)/(f_lower-f_upper) #updating the current value of xr
        f_xr=eq.subs(x,xr).evalf()
        a=f_lower*f_xr
        if i:
            err=abs(((xr-x0)/xr)*100) #Relative percentage error
            errors.append(err)
        else:
            errors.append(math.nan)
        if a==0:
            break
        elif a>0:
            lower=xr
            f_lower=f_xr
        else:
            upper=xr
            f_upper=f_xr
        if i:
            if (err<maxerror or f_xr==0):
                break
        x0=xr
    print(f"The root of the given funtion is {xr} ")
    plt.plot(errors,marker='o')
    plt.xlabel('Iteration number')
    plt.ylabel('Relative Error(%)')
    plt.title('Plot of relative errors vs the iteration number.')
    plt.grid()
    plt.show()
    plotFuntionClosed(equation,lo,up)
    return 0

def modified_false_position(equation):
    x=Symbol('x')
    eq=sympify(equation)
    lower=float(input('Please enter lower starting value (eg. 0.1)\n'))
    upper=float(input('Please enter upper starting value (eg. 1)\n'))
    lo=lower #declaring lo and up for use in 'plotFuntionClosed' funtion.
    up=upper
    f_lower=eq.subs(x,lower).evalf()
    f_upper=eq.subs(x,upper).evalf()
    if f_lower*f_upper>0:
        print("Error: Please select correct interval")
        quit()
    maxiters = int(input('Enter the Maximum Iterations:'))
    maxerror = float(input('Enter the Maximum Error (%):'))
    errors=[math.nan]
    il=0
    ul=0
    x0=lower
    for i in range(0,maxiters):
        xr=(f_lower*upper-f_upper*lower)/(f_lower-f_upper) #updating the current value of xr
        f_xr=eq.subs(x,xr).evalf()
        a=f_lower*f_xr
        if i:
            err=abs(((xr-x0)/xr)*100) #Relative percentage error
            errors.append(err)
        else:
            errors.append(math.nan)
        if a==0:
            break
        elif a>0:
            lower=xr
            f_lower=f_xr
            il=0
            ul=ul+1
            if ul>=2:
                f_upper=f_upper/2
        else:
            upper=xr
            f_upper=f_xr
            ul=0
            il=il+1
            if il>=2:
                f_lower=f_lower/2
        if i:
            if (err<maxerror or f_xr==0):
                break
        x0=xr
    print(f"The root of the given funtion is {xr} ")
    plt.plot(errors,marker='o')
    plt.xlabel('Iteration number')
    plt.ylabel('Relative Error(%)')
    plt.title('Plot of relative errors vs the iteration number.')
    plt.grid()
    plt.show()
    plotFuntionClosed(equation,lo,up)
    return 0

def newton_raphson(equation):
    x=Symbol('x')
    eq=sympify(equation)
    eqPrime=eq.diff(x,1)
    a=float(input('Please enter starting value (eg. 0.1)\n'))
    maxiters = int(input('Enter the Maximum Iterations:'))
    maxerror = float(input('Enter the Maximum Error (%):'))
    errors=[math.nan]
    b=a
    fx=eq.subs(x,a).evalf()
    for i in range(0,maxiters):
        b=a-fx/eqPrime.subs(x,a).evalf()
        if i:
            err=(abs(b-a)/b)*100 #Relative percentage error
            errors.append(err)
        else:
            errors.append(math.nan)
        fx=eq.subs(x,b).evalf()
        if i:
            if err<maxerror or fx==0 :
                break
        a=b
    print(f"The root of the given funtion is {a} ")
    plt.plot(errors,marker='o')
    plt.xlabel('Iteration number')
    plt.ylabel('Relative Error(%)')
    plt.title('Plot of relative errors vs the iteration number.')
    plt.grid()
    plt.show()
    plotFuntionOpen(equation)
    return 0

def secant(equation):
    x=Symbol('x')
    eq=sympify(equation)
    a=float(input('Please enter first starting value (eg. 0.1)\n'))
    b=float(input('Please enter  second starting value (eg. 0.1)\n'))
    fa=eq.subs(x,a).evalf()
    fb=eq.subs(x,b).evalf()
    maxiters = int(input('Enter the Maximum Iterations:'))
    maxerror = float(input('Enter the Maximum Error (%):'))
    errors=[math.nan]
    for i in range(0,maxiters):
        c=(fa*b-fb*a)/(fa-fb)
        fc=eq.subs(x,c).evalf()
        fa=fb
        fb=fc
        a=b
        b=c
        if i:
            err=(abs(a-b))/b*100 #Relative percentage error
            errors.append(err)
        if i:
            if err<maxerror or fb==0 :
                break
    plt.plot(errors,marker='o')
    plt.xlabel('Iteration number')
    plt.ylabel('Relative Error(%)')
    plt.title('Plot of relative errors vs the iteration number.')
    plt.grid()
    plt.show()
    print(b)
    plotFuntionOpen(equation)
    return 0

eq = input('Please enter a function (eg. 600*x**4 - 550*x**3 + 200*x**2 - 20*x - 1): \n')
method=input('Please select the method you want to use:\nBisection -> A\nFalse-position -> B\nModified false-position -> C\nNewton-Raphson -> D\nSecant -> E\n')

if method.upper()=='A':
    bisection(eq)
elif method.upper()=='B':
    false_position(eq)
elif method.upper()=='C':
    modified_false_position(eq)
elif method.upper()=='D':
    newton_raphson(eq)
elif method.upper()=='E':
    secant(eq)
else:
    print('Please select a correct option')