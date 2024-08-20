from IPython import get_ipython
get_ipython().run_line_magic('reset', '-sf')

import os
import shutil
import math
from sys import exit
from sympy import *
init_printing()
from sympy.solvers import solve
from mpmath import findroot

modelname = input('Enter the name of the system you would like to analyze: ')

if not os.path.isdir(modelname):
    print('The directory of that system was not found. Create it first and place ' + 
          'the data file inside named in the same way.')
    exit()
else:
    os.chdir(os.getcwd() + '\\' + modelname)

try:
    exec(open(modelname + '.py').read())
except:
    print('The file ' + modelname + '.py could not be run')
    exit()
    
nvar = len(var)
    
try:
    exec(open(os.path.dirname(os.path.realpath(__file__)) + '\\functions.py').read())
except:
    print('File functions.py is not in the same folder as the script you are running')
    exit()

file = open('Variables.txt','w')

try:
    var
except:
    print('Variables were not provided')
    exit()

for varnum in range(nvar):
    try:
        exec(var[varnum] + ' = symbols(var[varnum], real=True)')
        var[varnum]=eval(var[varnum])
    except:
        print('The script could not define your variable ' + (var[varnum]) + ' as a variable')
        exit()
    file.write(latex(var[varnum]) + '\n')
    
try:
    parameters
except:
    print('Parameters were not provided. The script will assume that there are no parameters.')
    parameters = []

npar=len(parameters)

if npar>0:
    newparameters=dict()
    for key in parameters.keys():
        try:
            exec(key + ' = symbols(key, real = True)')
            newparameters[eval(key)] = parameters[key]
        except:
            print('The script could not define your variable ' + key + ' as a variable')
            exit()
    parameters = newparameters
    
try:
    diffmatrix
except:
    print('The diffusion matrix was not provided.')
    exit()
    
try:
    for row in range(nvar):
        for col in range(nvar):
            diffmatrix[row][col] = eval(str(diffmatrix[row][col]).replace('^', '**'))
except:
    print('The diffusion matrix is not a function of the parameters of the system')
    exit()
        
diffmatrix = Matrix(diffmatrix)

# Kinetics

try:
    kinetics
except:
    print('Kinetics were not provided.')
    exit()
    
for functionnumber in range(nvar):
    try:
        kinetics[functionnumber] = eval(str(kinetics[functionnumber]).replace('^', '**'))
    except:
        print('The expression ' + kinetics[functionnumber] + ' is not a function of the parameters of your system')
        exit()

kinetics = Matrix(kinetics)

file = open('Kinetics.txt','w')
for functionnumber in range(nvar):
    file.write(latex(kinetics[functionnumber]) + '\n')
file.close()

jacobianmat = kinetics.jacobian(var)

# Equilibria

kineticsevaluated = kinetics
kineticsevaluated = kineticsevaluated.subs(parameters)
eq = solve(kineticsevaluated, var)

while eq==[]:
    print('This script could not find an equilibrium point for the given parameter values shown below: ')
    print(parameters)
    parameterchange = input('Would you like to change a parameter value? [y/n] ')
    while parameterchange!='y' and parameterchange!='n':
        parameterchange = input('You must enter your answer in the format [y/n]. ' +
                              'Would you like to change a parameter value? ')
    if parameterchange=='y':
        whichparam = input('Which parameter would you like to change? ')
        while True:
            try:
                whichparam = eval(whichparam)
                while whichparam not in parameters.keys():
                    whichparam = input('That is not a parameter of the system. ' +
                                     'Which parameter would you like to change? ')
                break
            except:
                whichparam = input('That is not a parameter of the system. Which parameter would you like to change? ')
        parameters[whichparam] = input('Enter a value of ' + whichparam + ': ')
        while not isfloat(parameters[whichparam]):
            parameters[whichparam] = input('What you entered before is not a number. ' +
                                         'Enter a value of ' + whichparam + ': ')
        parameters[whichparam] = eval(parameters[whichparam])
        kineticsevaluated = kinetics
        kineticsevaluated = kineticsevaluated.subs(parameters)
        eq = solve(kineticsevaluated, var)
    else:
        break
if isinstance(eq,list) and len(eq)>1:
    print('Your system has ' + str(len(eq)) + ' equilibria for the given parameter values given by:')
    for eqnum in range(len(eq)):
        print(str(eqnum+1) + '.- ' + str(eq[eqnum]) + '\n')
    eqnumber = input('Enter the number of the equilibrium point you want to consider: ')
    while not eqnumber.isnumeric() or int(eqnumber)==0:
        eqnumber = input('What you entered before is not a positive integer. ' +
                       'Enter the number of the equilibrium point you want to consider: ')
    eqnumber = int(eqnumber)
    eq=Matrix(list(eq[eqnumber-1]))        
elif isinstance(eq,list) and len(eq)==1:
    print('Your system has only one equilibrium point given by:')
    eqnumber = 0
    print(eq[0])
    eq = Matrix(eq[0])

# Wave conditions

muNF = symbols('mu_NF', real = True)
omegaNF = symbols('omega_NF', real = True)

dummyrealpart = symbols('dummyrealpart', real = True)

for counter1 in range(4):
    if counter1%2==0:
        second_set = range(0, 5, 2)
    else:
        second_set = range(1, 4, 2)
    for counter2 in second_set:
        if counter1==1 and counter2==1:
            coefmat11 = matrix('coefmat11', nvar, Add(jacobianmat, Mul(- 1, muNF, diffmatrix),
                                                      Mul(- 1, dummyrealpart, eye(nvar)),
                                                      Mul(- 1, I, omegaNF, eye(nvar))))
        else:
            exec(f'coefmat{counter1}{counter2} = matrix("coefmat{counter1}{counter2}", nvar, Add(jacobianmat, Mul(- Pow({counter1}, 2), muNF, diffmatrix), Mul(- {counter2}, I, omegaNF, eye(nvar))))')

# coefmat00 = matrix('coefmat00', nvar, jacobianmat)
# coefmat11 = matrix('coefmat11', nvar, Add(jacobianmat, Mul(- 1, muNF, diffmatrix),
#                                           Mul(- 1, dummyrealpart, eye(nvar)),
#                                           Mul(- 1, I, omegaNF, eye(nvar))))
# coefmat02 = matrix('coefmat02', nvar, Add(jacobianmat, Mul(- 2, I, omegaNF, eye(nvar))))
# coefmat20 = matrix('coefmat20', nvar, Add(jacobianmat, Mul(- 4, muNF, diffmatrix)))
# coefmat22 = matrix('coefmat22', nvar, Add(jacobianmat, Mul(- 4, muNF, diffmatrix),
#                                           Mul(-2, I, omegaNF, eye(nvar))))
# coefmat33 = matrix('coefmat33', nvar, Add(jacobianmat, Mul(- 9, muNF, diffmatrix),
#                                           Mul(- 3, I, omegaNF, eye(nvar))))
# coefmat31 = matrix('coefmat31', nvar, Add(jacobianmat, Mul(- 9, muNF, diffmatrix),
#                                           Mul(- 1, I, omegaNF, eye(nvar))))
# coefmat13 = matrix('coefmat13', nvar, Add(jacobianmat, Mul(- 1, muNF, diffmatrix),
#                                           Mul(- 3, I, omegaNF, eye(nvar))))
# coefmat04 = matrix('coefmat04', nvar, Add(jacobianmat, Mul(- 4, I, omegaNF, eye(nvar))))
# coefmat24 = matrix('coefmat24', nvar, Add(jacobianmat, Mul(- 4, muNF, diffmatrix),
#                                           Mul(- 4, I, omegaNF, eye(nvar))))

jacobianmatdet = coefmat11.dummy.det()

for row in range(nvar):
    for col in range(nvar):
        jacobianmatdet = jacobianmatdet.subs(coefmat11.dummy[row,col],coefmat11.actualcoord[row,col])

realdeterminant = expand(re(jacobianmatdet))
imagdeterminant = expand(im(jacobianmatdet))

coefmat11.actualcoord = coefmat11.actualcoord.subs(dummyrealpart, 0)

auxmat = matrix('dummy', 2)
auxmat.actualcoord = Matrix([realdeterminant, imagdeterminant]).jacobian([dummyrealpart,
                                                                    omegaNF]).subs(dummyrealpart, 0)

auxiliaryterm = auxmat.dummy.det()

for row in range(2):
    for col in range(2):
        auxiliaryterm = auxiliaryterm.subs(auxmat.dummy[row, col], auxmat.actualcoord[row, col])

realdeterminant = realdeterminant.subs(dummyrealpart, 0)
imagdeterminant = imagdeterminant.subs(dummyrealpart, 0)

auxmat.actualcoord = Matrix([realdeterminant, imagdeterminant]).jacobian([muNF, omegaNF])

realdeterminantderivative = auxmat.dummy.det()

for row in range(2):
    for col in range(2):
        realdeterminantderivative = realdeterminantderivative.subs(auxmat.dummy[row, col], auxmat.actualcoord[row, col])
        
realdeterminanteval = realdeterminant
imagdeterminanteval = imagdeterminant
derivativeeval = realdeterminantderivative

# if eq!=[]:
#     for varnum in range(nvar):
#         try:
#             determinanteval = determinanteval.subs(var[varnum], eq[varnum])
#             derivativeeval = derivativeeval.subs(var[varnum], eq[varnum])
#         except:
#             determinanteval = determinanteval.subs(var[varnum], eq[var[varnum]])
#             derivativeeval = derivativeeval.subs(var[varnum], eq[var[varnum]])
#     determinanteval = determinanteval.subs(parameters)
#     derivativeeval = derivativeeval.subs(parameters)
#     tol = 7
#     imtol = 15
#     getout = 0
#     try:
#         mucritical = solve(derivativeeval,muNF)
#         for muvalue in mucritical:
#             if (abs(simplify(determinanteval.subs(muNF,muvalue)))<5*10**(-tol)
#                 and abs(complex(muvalue).imag)<10**(-imtol)):
#                 ksquared = complex(muvalue).real
#                 getout = 1
#                 break
#     except:
#         pass
#     if getout==0:
#         determinanteval = jacobianmatdet
#         derivativeeval = determinantderivative
#         parametertochange = input('There is no Wave bifurcation for the parameters provided. Enter a parameter that can be changed in order to find the wavenumber: ')
#         counter = 0
#         while True:
#             try:
#                 if counter==1:
#                     parametertochange = input('Enter a parameter that can be changed in order to find the wavenumber: ')
#                 parametertochange = eval(parametertochange)
#                 if parametertochange not in parameters.keys():
#                     print('You did not enter a parameter of the system.')
#                     continue
#                 else:
#                     if parametertochange in parameters.keys():
#                         if (determinanteval==determinanteval.subs(parametertochange,parameters[parametertochange])
#                             and derivativeeval==derivativeeval.subs(parametertochange,parameters[parametertochange])):
#                             print('Your parameter does not produce any changes in the equations for the Wave bifurcation.')
#                             continue
#                 break
#             except:
#                 counter = 1
#         kineticsevaluated=kinetics
#         for key in parameters:
#             if key!=parametertochange:
#                 kineticsevaluated = kineticsevaluated.subs(key,parameters[key])
#                 determinanteval = determinanteval.subs(key,parameters[key])
#                 derivativeeval = derivativeeval.subs(key,parameters[key])
#             else:
#                 initialpartochange = parameters[key]
#         initialeq = eq
#         eq = solve(kineticsevaluated,var)[eqnumber]
#         for varnum in range(nvar):
#             determinanteval = determinanteval.subs(var[varnum],eq[varnum])
#             derivativeeval = derivativeeval.subs(var[varnum],eq[varnum])
#         try:
#             zerofunction = [lambda mutofind, parametertofind: determinanteval.subs([(muNF,mutofind),(parametertochange,parametertofind)]), lambda mutofind, parametertofind: derivativeeval.subs([(muNF,mutofind),(parametertochange,parametertofind)])]
#             [muvalue,newparval] = findroot(zerofunction, (1,initialpartochange))
#             if abs(simplify(determinanteval.subs([(muNF,muvalue),(parametertochange,newparval)])))<10**(-tol) and abs(simplify(derivativeeval.subs([(muNF,muvalue),(parametertochange,newparval)])))<10**(-tol) and abs(complex(muvalue).imag)<10**(-imtol) and abs(complex(newparval).imag)<10**(-imtol):
#                 ksquared=complex(muvalue).real
#                 parameters[parametertochange]=complex(newparval).real
#                 getout=1
#         except:
#             pass
#         if getout==0:
#             try:
#                 intersection=resultant(determinanteval,derivativeeval,muNF)
#                 newparval=float(findroot(lambda parametertofind: intersection.subs(parametertochange,parametertofind),initialpartochange))
#                 mucritical=solve(derivativeeval.subs(parametertochange,newparval),muNF)
#                 for muvalue in mucritical:
#                     if abs(simplify(determinanteval.subs([(muNF,muvalue),(parametertochange,newparval)])))<10**(-tol) and abs(complex(muvalue).imag)<10**(-imtol) and abs(complex(newparval).imag)<10**(-imtol):
#                         ksquared=complex(muvalue).real
#                         getout = 1
#                         parameters[parametertochange] = complex(newparval).real
#                         break
#             except:
#                 pass
#             if getout==0:
#                 try:
#                     mucritical = solve(derivativeeval,muNF)
#                     for muvalue in mucritical:
#                         try:
#                             newparval = float(findroot(lambda parametertofind: determinanteval.subs([(muNF,muvalue),(parametertochange,parametertofind)]),initialpartochange))
#                             if simplify(determinanteval.subs([(muNF,muvalue),(parametertochange,newparval)]))<10**(-tol) and complex(simplify(muvalue.subs(parametertochange,newparval))).imag<10**(-imtol) and complex(simplify(newparval)).imag<10**(-imtol):
#                                 ksquared = complex(simplify(muvalue.subs(parametertochange,newparval))).real
#                                 getout = 1
#                                 parameters[parametertochange] = complex(newparval).real
#                                 break
#                         except:
#                             continue
#                 except:
#                     pass
#                 if getout==0:
#                     eq = initialeq
#                     kval = input('This script could not find the value of k. Make sure that you have a Wave bifurcation for the parameter values provided. If you are at a Wave bifurcation. Enter a value of k you want to consider: ')
#                     while not isfloat(kval):
#                         kval = input('What you entered before is not a number. Enter a value of k you want to consider: ')
#                     kval = eval(kval)
#                     ksquared = Pow(kval, 2)
# else:
#     ksquared = 0

tol = 1e-7
ksquared = 0
omegaval = 0
    
file = open('Real part of the determinant of the Jacobian matrix.txt','w')
file.write(latex(realdeterminant))
file.close()
file = open('Imaginary part of the determinant of the Jacobian matrix.txt','w')
file.write(latex(imagdeterminant))
file.close()
file = open('Real part of the derivative of the Determinant.txt','w')
file.write(latex(realdeterminantderivative))
file.close()
file = open("Non-zero variable.txt",'w')
file.write(latex(auxiliaryterm))
file.close()

# Normal form

negativeRHS = Vector('dummynegativeRHS')

firstorderderivatives = list()
secondorderderivatives = list()
thirdorderderivatives = list()
fourthorderderivatives = list()
fifthorderderivatives = list()
for counter1 in range(nvar):
    firstorderderivatives.append(diff(kinetics, var[counter1]))
    secondorderderivatives.append(list())
    thirdorderderivatives.append(list())
    fourthorderderivatives.append(list())
    fifthorderderivatives.append(list())
    for counter2 in range(nvar):
        secondorderderivatives[counter1].append(diff(firstorderderivatives[counter1], var[counter2]))
        thirdorderderivatives[counter1].append(list())
        fourthorderderivatives[counter1].append(list())
        fifthorderderivatives[counter1].append(list())
        for counter3 in range(nvar):
            thirdorderderivatives[counter1][counter2].append(diff(secondorderderivatives[counter1][counter2],
                                                                  var[counter3]))
            fourthorderderivatives[counter1][counter2].append(list())
            fifthorderderivatives[counter1][counter2].append(list())
            for counter4 in range(nvar):
                fourthorderderivatives[counter1][counter2][counter3].append(
                    diff(thirdorderderivatives[counter1][counter2][counter3], var[counter4]))
                fifthorderderivatives[counter1][counter2][counter3].append(list())
                for counter5 in range(nvar):
                    fifthorderderivatives[counter1][counter2][counter3][counter4].append(
                        diff(fourthorderderivatives[counter1][counter2][counter3][counter4], var[counter5]))

phi11NF = Vector('phi11NF')
psi11NF = Vector('psi11NF')
W02NF = Vector('W02NF')
W22NF = Vector('W22NF')
W011NF = Vector('W011NF')
W211NF = Vector('W211NF')
W13NF = Vector('W13NF')
W33NF = Vector('W33NF')
W121NF = Vector('W121NF')
W1221NF = Vector('W1221NF')
W321NF = Vector('W321NF')
W04NF = Vector('W04NF')
W24NF = Vector('W24NF')
W031NF = Vector('W031NF')
W231NF = Vector('W231NF')
W022NF = Vector('W022NF')
W222NF = Vector('W222NF')

getout = 0

for row in range(nvar):
    for col in range(nvar):
        submatrixrows = list(range(nvar))
        submatrixcols = list(range(nvar))
        submatrixrows.remove(row)
        submatrixcols.remove(col)
        invertiblesubmatrix = coefmat11.actualcoord.extract(submatrixrows, submatrixcols)
        submatrixeval = invertiblesubmatrix
        submatrixeval = submatrixeval.subs(parameters)
        if eq!=[]:
            for varnum in range(nvar):
                submatrixeval = submatrixeval.subs(var[varnum], eq[varnum])
            submatrixeval = submatrixeval.subs(muNF, ksquared)
            submatrixeval = submatrixeval.subs(omegaNF, omegaval)
            if abs(N(submatrixeval.det())) > tol:
                criticalrow = row
                criticalcol = col
                getout = 1
                break
        else:
            criticalrow = 0
            criticalcol = 0
            submatrixrows = list(range(nvar))
            submatrixcols = list(range(nvar))
            submatrixrows.remove(0)
            submatrixcols.remove(0)
            invertiblesubmatrix = coefmat11.actualcoord.extract(submatrixrows, submatrixcols)
            break
    if getout==1:
        break
    
coefsubmatrix = matrix('dummysubmatrix', nvar-1, invertiblesubmatrix)

phi11NF = kernel_determination(phi11NF, coefmat11, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

if phiunit=='y':
    phi11NF.actualcoord = Mul(Pow(sqrt(phi11NF.actualcoord.dot(conjugate(phi11NF.actualcoord))),-1),
                              phi11NF.actualcoord)
    
phi11NF_eval = evaluation_dict(phi11NF)

print('First order ready')

DS_phi11conjphi11 = secondorderapplied(phi11NF, phi11NF.conj)
DS_phi11phi11 = secondorderapplied(phi11NF, phi11NF)

negativeRHS.actualcoord = DS_phi11conjphi11

W02NF = linearsolver(W02NF, negativeRHS, coefmat00)

negativeRHS.actualcoord = DS_phi11phi11

W22NF = linearsolver(W22NF, negativeRHS, coefmat22)

negativeRHS.actualcoord = Mul(2, DS_phi11phi11)

W011NF = linearsolver(W011NF, negativeRHS, coefmat02)

negativeRHS.actualcoord = Mul(2, DS_phi11conjphi11)

W211NF = linearsolver(W211NF, negativeRHS, coefmat20)

W02NF_eval = evaluation_dict(W02NF)
W22NF_eval = evaluation_dict(W22NF)
W011NF_eval = evaluation_dict(W011NF)
W211NF_eval = evaluation_dict(W211NF)
        
print('Second order ready')

DS_phi11W02 = secondorderapplied(phi11NF, W02NF)
DS_conjphi11W22 = secondorderapplied(phi11NF.conj, W22NF)

DS_conjphi11W011 = secondorderapplied(phi11NF.conj, W011NF)
DS_phi11W211 = secondorderapplied(phi11NF, W211NF)

TS_phi11phi11conjphi11 = thirdorderapplied(phi11NF, phi11NF, phi11NF.conj)

Tcoefsubmatrix = matrix('Tdummysubmatrix', nvar - 1, transpose(invertiblesubmatrix))
Tcoefmat11 = matrix('Tcoefmat', nvar, transpose(coefmat11.actualcoord))

psi11NF = kernel_determination(psi11NF, Tcoefmat11, criticalrow, Tcoefsubmatrix, submatrixcols, submatrixrows)

psi11NF_eval = evaluation_dict(psi11NF)

denominator = phi11NF.dummy.dot(psi11NF.dummy)

C130 = Mul(Pow(denominator, -1), psi11NF.dummy.dot(Add(Mul(4, DS_phi11W02), Mul(2, DS_conjphi11W22),
                                                       Mul(3, TS_phi11phi11conjphi11))))

C112 = Mul(Pow(denominator, -1),
           psi11NF.dummy.dot(Add(Mul(2, DS_conjphi11W011), Mul(4, DS_phi11W02), Mul(2, DS_phi11W211),
                                Mul(6, TS_phi11phi11conjphi11))))
    
C130conj = conjugate(C130)

C112conj = conjugate(C112)

if crosscoef=='y':
    try:    
        crosspar = eval(crosspar)
        
        pre_C11 = Mul(Pow(denominator, -1),
                      psi11NF.dummy.dot(Mul(Add(jacobianmat, Mul(- muNF, diffmatrix)), phi11NF.dummy)))
        
        # file = open('Kernels.txt', 'w')
        # for varnum in range(nvar):
        #     if varnum<nvar - 1:
        #         file.write(latex(phi11NF.actualcoord[varnum]) + ',')
        #     else:
        #         file.write(latex(phi11NF.actualcoord[varnum]) + '\n')
        # for varnum in range(nvar):
        #     if varnum<nvar - 1:
        #         file.write(latex(psi11NF.actualcoord[varnum]) + ',')
        #     else:
        #         file.write(latex(psi11NF.actualcoord[varnum]))
        # file.close()

        file = open('To get cross-order coefficient.txt', 'w')
        file.write(latex(crosspar) + '\n')
        file.write(latex(pre_C11))
        file.close()
    except:
        print('The parameter crosspar is not a parameter of the system.')

if fifthcoef=='y':
    DS_phi11W22 = secondorderapplied(phi11NF, W22NF)
    DS_phi11W011 = secondorderapplied(phi11NF, W011NF)
    
    TS_phi11phi11phi11 = thirdorderapplied(phi11NF, phi11NF, phi11NF)
    
    negativeRHS.actualcoord = Add(Mul(-1, C130, phi11NF.dummy), Mul(4, DS_phi11W02), Mul(2, DS_conjphi11W22),
                                  Mul(3, TS_phi11phi11conjphi11))
    
    W13NF = critical_linearsolver(W13NF, negativeRHS, criticalcol, coefsubmatrix,
                                  submatrixrows, submatrixcols)

    if orthogonal=='y':
        if phiunit=='y':
            W13NF.actualcoord = Add(W13NF.actualcoord, Mul(-1, W13NF.actualcoord.dot(phi11NF.conj.dummy),
                                                            phi11NF.dummy))
        else:
            W13NF.actualcoord = Add(W13NF.actualcoord, Mul(-1, W13NF.actualcoord.dot(phi11NF.conj.dummy),
                                                            Pow(phi11NF.dummy.dot(phi11NF.conj.dummy), -1),
                                                            phi11NF.dummy))
    
    negativeRHS.actualcoord = Add(Mul(-1, C112, phi11NF.dummy), Mul(2, DS_conjphi11W011), Mul(4, DS_phi11W02),
                                  Mul(2, DS_phi11W211), Mul(6, TS_phi11phi11conjphi11))
    
    W121NF = critical_linearsolver(W121NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)
                
    if orthogonal=='y':
        if phiunit=='y':
            W121NF.actualcoord = Add(W121NF.actualcoord, Mul(-1, W121NF.actualcoord.dot(phi11NF.conj.dummy),
                                                            phi11NF.dummy))
        else:
            W121NF.actualcoord = Add(W121NF.actualcoord, Mul(-1, W121NF.actualcoord.dot(phi11NF.conj.dummy),
                                                            Pow(phi11NF.dummy.dot(phi11NF.conj.dummy), -1),
                                                            phi11NF.dummy))
        
    negativeRHS.actualcoord = Add(Mul(2, DS_phi11W22), TS_phi11phi11phi11)
    
    W33NF = linearsolver(W33NF, negativeRHS, coefmat33)
    
    negativeRHS.actualcoord = Add(Mul(2, DS_phi11W22), Mul(2, DS_phi11W011), Mul(3, TS_phi11phi11phi11))
    
    W1221NF = linearsolver(W1221NF, negativeRHS, coefmat13)
    
    negativeRHS.actualcoord = Add(Mul(2, DS_phi11W211), Mul(2, DS_conjphi11W22), Mul(3, TS_phi11phi11conjphi11))
    
    W321NF = linearsolver(W321NF, negativeRHS, coefmat31)
    
    W13NF_eval = evaluation_dict(W13NF)
    W121NF_eval = evaluation_dict(W121NF)
    W33NF_eval = evaluation_dict(W33NF)
    W1221NF_eval = evaluation_dict(W1221NF)
    W321NF_eval = evaluation_dict(W321NF)
    
    print('Third order ready')
    
    DS_W02W02 = secondorderapplied(W02NF, W02NF)
    DS_W22conjW22 = secondorderapplied(W22NF, W22NF.conj)
    DS_phi11conjW13 = secondorderapplied(phi11NF, W13NF.conj)
    
    DS_W02W22 = secondorderapplied(W02NF, W22NF)
    DS_phi11W13 = secondorderapplied(phi11NF, W13NF)
    DS_conjphi11W33 = secondorderapplied(phi11NF.conj, W33NF)
    
    DS_W02W011 = secondorderapplied(W02NF, W011NF)
    DS_W22W211 = secondorderapplied(W22NF, W211NF)
    DS_phi11W121 = secondorderapplied(phi11NF, W121NF)
    DS_conjphi11W1221 = secondorderapplied(phi11NF.conj, W1221NF)
    
    DS_W02W211 = secondorderapplied(W02NF, W211NF)
    DS_W22conjW011 = secondorderapplied(W22NF, W011NF.conj)
    DS_phi11conjW121 = secondorderapplied(phi11NF, W121NF.conj)
    DS_conjphi11W13 = secondorderapplied(phi11NF.conj, W13NF)
    DS_conjphi11W321 = secondorderapplied(phi11NF.conj, W321NF)
    
    DS_W211W211 = secondorderapplied(W211NF, W211NF)
    DS_W011conjW011 = secondorderapplied(W011NF, W011NF.conj)
    
    DS_W011W211 = secondorderapplied(W011NF, W211NF)
    DS_phi11W321 = secondorderapplied(phi11NF, W321NF)
    
    TS_phi11phi11conjW22 = thirdorderapplied(phi11NF, phi11NF, W22NF.conj)
    TS_phi11conjphi11W02 = thirdorderapplied(phi11NF, phi11NF.conj, W02NF)
    
    TS_phi11phi11W02 = thirdorderapplied(phi11NF, phi11NF, W02NF)
    TS_phi11conjphi11W22 = thirdorderapplied(phi11NF, phi11NF.conj, W22NF)
    
    TS_phi11phi11W211 = thirdorderapplied(phi11NF, phi11NF, W211NF)
    TS_phi11conjphi11W011 = thirdorderapplied(phi11NF, phi11NF.conj, W011NF)
        
    TS_phi11phi11conjW011 = thirdorderapplied(phi11NF, phi11NF, W011NF.conj)
    TS_phi11conjphi11W211 = thirdorderapplied(phi11NF, phi11NF.conj, W211NF)
    TS_conjphi11conjphi11W22 = thirdorderapplied(phi11NF.conj, phi11NF.conj, W22NF)
    
    Q4S_phi11phi11conjphi11conjphi11 = fourthorderapplied(phi11NF, phi11NF, phi11NF.conj, phi11NF.conj)
    Q4S_phi11phi11phi11conjphi11 = fourthorderapplied(phi11NF, phi11NF, phi11NF, phi11NF.conj)
    
    negativeRHS.actualcoord = Add(Mul(-2, C130, W02NF.dummy), Mul(2, DS_W02W02),
                                  DS_W22conjW22, Mul(2, DS_phi11conjW13),
                                  Mul(3, TS_phi11phi11conjW22), Mul(6, TS_phi11conjphi11W02),
                                  Mul(3, Q4S_phi11phi11conjphi11conjphi11))
    
    W04NF = linearsolver(W04NF, negativeRHS, coefmat00)
    
    negativeRHS.actualcoord = Add(Mul(-2, C130, W22NF.dummy), Mul(4, DS_W02W22), Mul(2, DS_phi11W13),
                                  Mul(2, DS_conjphi11W33), Mul(6, TS_phi11phi11W02), Mul(6, TS_phi11conjphi11W22),
                                  Mul(4, Q4S_phi11phi11phi11conjphi11))
    
    W24NF = linearsolver(W24NF, negativeRHS, coefmat22)
    
    negativeRHS.actualcoord = Add(Mul(-1, Add(C130, C112), W011NF.dummy), Mul(4, DS_W02W011), Mul(2, DS_W22W211),
                                  Mul(2, DS_phi11W13), Mul(2, DS_phi11W121), Mul(2, DS_conjphi11W1221),
                                  Mul(12, TS_phi11phi11W02), Mul(3, TS_phi11phi11W211),
                                  Mul(6, TS_phi11conjphi11W22), Mul(6, TS_phi11conjphi11W011),
                                  Mul(12, Q4S_phi11phi11phi11conjphi11))
    
    W031NF = linearsolver(W031NF, negativeRHS, coefmat02)
    
    negativeRHS.actualcoord = Add(Mul(-1, Add(C130, C112conj), W211NF.dummy),
                                  Mul(4, DS_W02W211), Mul(2, DS_W22conjW011), Mul(2, DS_phi11conjW121),
                                  Mul(2, DS_conjphi11W13), Mul(2, DS_conjphi11W321),
                                  Mul(3, TS_phi11phi11conjW011), Mul(12, TS_phi11conjphi11W02),
                                  Mul(6, TS_phi11conjphi11W211), Mul(6, TS_conjphi11conjphi11W22),
                                  Mul(12, Q4S_phi11phi11conjphi11conjphi11))
    
    W231NF = linearsolver(W231NF, negativeRHS, coefmat20)
    
    negativeRHS.actualcoord = Add(Mul(-4, C112, W02NF.dummy), Mul(4, DS_W02W02),
                                  DS_W211W211, DS_W011conjW011, Mul(4, DS_phi11conjW121),
                                  Mul(12, TS_phi11conjphi11W02), Mul(6, TS_phi11conjphi11W211),
                                  Mul(6, TS_phi11phi11conjW011), Mul(12, Q4S_phi11phi11conjphi11conjphi11))
        
    W022NF = linearsolver(W022NF, negativeRHS, coefmat00)
    
    negativeRHS.actualcoord = Add(Mul(-2, C112, W22NF.dummy), Mul(2, DS_W011W211), Mul(4, DS_W02W22),
                                  Mul(2, DS_phi11W121), Mul(2, DS_phi11W321),
                                  Mul(2, DS_conjphi11W1221), Mul(6, TS_phi11phi11W02),
                                  Mul(6, TS_phi11phi11W211), Mul(6, TS_phi11conjphi11W22),
                                  Mul(6, TS_phi11conjphi11W011), Mul(12, Q4S_phi11phi11phi11conjphi11))
    
    W222NF = linearsolver(W222NF, negativeRHS, coefmat22)
    
    W04NF_eval = evaluation_dict(W04NF)
    W24NF_eval = evaluation_dict(W24NF)
    W031NF_eval = evaluation_dict(W031NF)
    W231NF_eval = evaluation_dict(W231NF)
    W022NF_eval = evaluation_dict(W022NF)
    W222NF_eval = evaluation_dict(W222NF)
    
    print('Fourth order ready')
    
    DS_phi11W04 = secondorderapplied(phi11NF, W04NF)
    DS_phi11conjW04 = secondorderapplied(phi11NF, W04NF.conj)
    DS_conjphi11W24 = secondorderapplied(phi11NF.conj, W24NF)
    DS_W02W13 = secondorderapplied(W02NF, W13NF)
    DS_W22conjW13 = secondorderapplied(W22NF, W13NF.conj)
    DS_conjW22W33 = secondorderapplied(W22NF.conj, W33NF)
    
    DS_phi11conjW231 = secondorderapplied(phi11NF, W231NF.conj)
    DS_conjphi11W031 = secondorderapplied(phi11NF.conj, W031NF)
    DS_W02W121 = secondorderapplied(W02NF, W121NF)
    DS_conjW22W1221 = secondorderapplied(W22NF.conj, W1221NF)
    DS_W22conjW321 = secondorderapplied(W22NF, W321NF.conj)
    DS_W13W211 = secondorderapplied(W13NF, W211NF)
    DS_conjW13W011 = secondorderapplied(W13NF.conj, W011NF)
    
    DS_phi11W022 = secondorderapplied(phi11NF, W022NF)
    DS_phi11conjW022 = secondorderapplied(phi11NF, W022NF.conj)
    DS_phi11W231 = secondorderapplied(phi11NF, W231NF)
    DS_conjphi11W222 = secondorderapplied(phi11NF.conj, W222NF)
    DS_conjW121W22 = secondorderapplied(W121NF.conj, W22NF)
    DS_conjW121W011 = secondorderapplied(W121NF.conj, W011NF)
    DS_W211W121 = secondorderapplied(W211NF, W121NF)
    DS_W211W321 = secondorderapplied(W211NF, W321NF)
    DS_conjW011W1221 = secondorderapplied(W011NF.conj, W1221NF)
    
    TS_phi11W02W02 = thirdorderapplied(phi11NF, W02NF, W02NF)
    TS_phi11W22conjW22 = thirdorderapplied(phi11NF, W22NF, W22NF.conj)
    TS_conjphi11W02W22 = thirdorderapplied(phi11NF.conj, W02NF, W22NF)
    TS_phi11phi11conjW13 = thirdorderapplied(phi11NF, phi11NF, W13NF.conj)
    TS_phi11conjphi11W13 = thirdorderapplied(phi11NF, phi11NF.conj, W13NF)
    TS_conjphi11conjphi11W33 = thirdorderapplied(phi11NF.conj, phi11NF.conj, W33NF)
    
    TS_phi11W02W211 = thirdorderapplied(phi11NF, W02NF, W211NF)
    TS_phi11conjW22W22 = thirdorderapplied(phi11NF, W22NF.conj, W22NF)
    TS_phi11conjW22W011 = thirdorderapplied(phi11NF, W22NF.conj, W011NF)
    TS_conjphi11W02W011 = thirdorderapplied(phi11NF.conj, W02NF, W011NF)
    TS_conjphi11W22W211 = thirdorderapplied(phi11NF.conj, W22NF, W211NF)
    TS_phi11conjphi11W121 = thirdorderapplied(phi11NF, phi11NF.conj, W121NF)
    TS_conjphi11conjphi11W1221 = thirdorderapplied(phi11NF.conj, phi11NF.conj, W1221NF)
    TS_phi11phi11conjW321 = thirdorderapplied(phi11NF, phi11NF, W321NF.conj)
    
    TS_phi11W211W211 = thirdorderapplied(phi11NF, W211NF, W211NF)
    TS_phi11conjW011W22 = thirdorderapplied(phi11NF, W011NF.conj, W22NF)
    TS_phi11conjW011W011 = thirdorderapplied(phi11NF, W011NF.conj, W011NF)
    TS_conjphi11W211W22 = thirdorderapplied(phi11NF.conj, W211NF, W22NF)
    TS_conjphi11W211W011 = thirdorderapplied(phi11NF.conj, W211NF, W011NF)
    TS_phi11phi11conjW121 = thirdorderapplied(phi11NF, phi11NF, W121NF.conj)    
    TS_phi11conjphi11W321 = thirdorderapplied(phi11NF, phi11NF.conj, W321NF)
    
    Q4S_phi11phi11phi11conjW22 = fourthorderapplied(phi11NF, phi11NF, phi11NF, W22NF.conj)
    Q4S_phi11phi11conjphi11W02 = fourthorderapplied(phi11NF, phi11NF, phi11NF.conj, W02NF)
    Q4S_phi11conjphi11conjphi11W22 = fourthorderapplied(phi11NF, phi11NF.conj, phi11NF.conj, W22NF)
    
    Q4S_phi11phi11conjphi11W211 = fourthorderapplied(phi11NF, phi11NF, phi11NF.conj, W211NF)
    Q4S_phi11conjphi11conjphi11W011 = fourthorderapplied(phi11NF, phi11NF.conj, phi11NF.conj, W011NF)
    
    Q4S_phi11phi11phi11conjW011 = fourthorderapplied(phi11NF, phi11NF, phi11NF, W011NF.conj)
    
    Q5S_phi11phi11phi11conjphi11conjphi11 = fifthorderapplied(phi11NF, phi11NF, phi11NF, phi11NF.conj, phi11NF.conj)
        
    C150 = Mul(Pow(denominator, -1),
                psi11NF.dummy.dot(Add(Mul(-1, Add(Mul(2, C130), C130conj), W13NF.dummy),
                                      Mul(2, DS_phi11W04), Mul(2, DS_phi11conjW04), Mul(2, DS_conjphi11W24),
                                      Mul(4, DS_W02W13), Mul(2, DS_W22conjW13), Mul(2, DS_conjW22W33),
                                      Mul(12, TS_phi11W02W02), Mul(6, TS_phi11W22conjW22),
                                      Mul(12, TS_conjphi11W02W22), Mul(3, TS_phi11phi11conjW13),
                                      Mul(6, TS_phi11conjphi11W13), Mul(3, TS_conjphi11conjphi11W33),
                                      Mul(4, Q4S_phi11phi11phi11conjW22), Mul(24, Q4S_phi11phi11conjphi11W02),
                                      Mul(12, Q4S_phi11conjphi11conjphi11W22),
                                      Mul(10, Q5S_phi11phi11phi11conjphi11conjphi11))))
    
    C114 = Mul(Pow(denominator, -1),
                psi11NF.dummy.dot(Add(Mul(-1, Add(C112, C130, C130conj), W121NF.dummy),
                                      Mul(2, DS_phi11W04), Mul(2, DS_phi11conjW04), Mul(2, DS_phi11conjW231),
                                      Mul(2, DS_conjphi11W031), Mul(4, DS_W02W121), Mul(2, DS_conjW22W1221),
                                      Mul(2, DS_W22conjW321), Mul(2, DS_W13W211), Mul(2, DS_conjW13W011),
                                      Mul(12, Add(TS_phi11W02W02, TS_phi11W02W211)), Mul(6, TS_phi11conjW22W22),
                                      Mul(6, TS_phi11conjW22W011), Mul(12, TS_conjphi11W02W011),
                                      Mul(6, TS_conjphi11W22W211), Mul(6, Add(TS_phi11conjphi11W13,
                                      TS_phi11conjphi11W121)), Mul(3, TS_conjphi11conjphi11W1221),
                                      Mul(6, TS_phi11phi11conjW13), Mul(3, TS_phi11phi11conjW321),
                                      Mul(48, Q4S_phi11phi11conjphi11W02), Mul(12, Q4S_phi11phi11conjphi11W211),
                                      Mul(12, Q4S_phi11phi11phi11conjW22),
                                      Mul(12, Add(Q4S_phi11conjphi11conjphi11W22, Q4S_phi11conjphi11conjphi11W011)),
                                      Mul(30, Q5S_phi11phi11phi11conjphi11conjphi11))))
        
    C132 = Mul(Pow(denominator, -1),
                psi11NF.dummy.dot(Add(Mul(-1, Add(C130, C112, C112conj), W121NF.dummy),
                                      Mul(-1, Add(Mul(2, C112), C112conj), W13NF.dummy), Mul(2, DS_phi11W022),
                                      Mul(2, DS_phi11conjW022), Mul(2, DS_phi11W231), Mul(2, DS_conjphi11W222),
                                      Mul(2, DS_conjphi11W031), Mul(4, DS_W02W121), Mul(2, DS_conjW121W22),
                                      Mul(2, DS_conjW121W011), Mul(2, DS_W211W121), Mul(2, DS_W211W321),
                                      Mul(2, DS_conjW011W1221), Mul(4, DS_W02W13), Mul(6, TS_phi11W211W211),
                                      Mul(6, TS_phi11conjW011W22), Mul(6, TS_phi11conjW011W011),
                                      Mul(12, Add(TS_conjphi11W02W22, TS_conjphi11W02W011)),
                                      Mul(6, Add(TS_conjphi11W211W22, TS_conjphi11W211W011)),
                                      Mul(12, Add(Mul(2, TS_phi11W02W02), TS_phi11W02W211)),
                                      Mul(9, TS_phi11phi11conjW121), Mul(6, TS_phi11conjphi11W13),
                                      Mul(12, TS_phi11conjphi11W121), Mul(6, TS_phi11conjphi11W321),
                                      Mul(6, TS_conjphi11conjphi11W1221),
                                      Mul(36, Add(Mul(2, Q4S_phi11phi11conjphi11W02),
                                                  Q4S_phi11phi11conjphi11W211)),
                                      Mul(12, Q4S_phi11phi11phi11conjW011),
                                      Mul(24, Q4S_phi11conjphi11conjphi11W22),
                                      Mul(24, Q4S_phi11conjphi11conjphi11W011),
                                      Mul(60, Q5S_phi11phi11phi11conjphi11conjphi11))))
    
    print('The calculation of the fifth-order coefficients was carried out successfully.')
    
# C130 = C130.subs(W02NF_eval)
# C130 = C130.subs(W22NF_eval)
# C130 = C130.subs(W011NF_eval)
# C130 = C130.subs(W211NF_eval)

# C112 = C112.subs(W02NF_eval)
# C112 = C112.subs(W22NF_eval)
# C112 = C112.subs(W011NF_eval)
# C112 = C112.subs(W211NF_eval)

# C130 = C130.subs(phi11NF_eval)
# C130 = C130.subs(psi11NF_eval)

# C112 = C112.subs(phi11NF_eval)
# C112 = C112.subs(psi11NF_eval)

C130 = simplify(C130)
C112 = simplify(C112)
    
file = open('Third-order coefficient C130.txt', 'w')
file.write(latex(C130))
file.close()

file = open('Third-order coefficient C112.txt', 'w')
file.write(latex(C112))
file.close()

file = open('Vectors.txt', 'w')
write_vector(phi11NF, file)
write_vector(psi11NF, file)
write_vector(W02NF, file)
write_vector(W22NF, file)
write_vector(W011NF, file)
write_vector(W211NF, file)
file.close()

print('The third-order coefficients were computed and saved into text files.')
    
if fifthcoef=='n':
    print('Everything but the fifth-order coefficients was computed and saved into files.')
elif fifthcoef=='y':
    file = open('Vectors.txt', 'a')
    write_vector(W13NF, file)
    write_vector(W33NF, file)
    write_vector(W121NF, file)
    write_vector(W1221NF, file)
    write_vector(W321NF, file)
    
    write_vector(W04NF, file)
    write_vector(W24NF, file)
    write_vector(W031NF, file)
    write_vector(W231NF, file)
    write_vector(W022NF, file)
    write_vector(W222NF, file)
    
    file.close()
    
    # C150 = C150.subs(W04NF_eval)
    # C150 = C150.subs(W24NF_eval)
    # C150 = C150.subs(W031NF_eval)
    # C150 = C150.subs(W231NF_eval)
    # C150 = C150.subs(W022NF_eval)
    # C150 = C150.subs(W222NF_eval)
    
    # C114 = C114.subs(W04NF_eval)
    # C114 = C114.subs(W24NF_eval)
    # C114 = C114.subs(W031NF_eval)
    # C114 = C114.subs(W231NF_eval)
    # C114 = C114.subs(W022NF_eval)
    # C114 = C114.subs(W222NF_eval)
    
    # C132 = C132.subs(W04NF_eval)
    # C132 = C132.subs(W24NF_eval)
    # C132 = C132.subs(W031NF_eval)
    # C132 = C132.subs(W231NF_eval)
    # C132 = C132.subs(W022NF_eval)
    # C132 = C132.subs(W222NF_eval)
    
    # C150 = C150.subs(W13NF_eval)
    # C150 = C150.subs(W33NF_eval)
    # C150 = C150.subs(W121NF_eval)
    # C150 = C150.subs(W1221NF_eval)
    # C150 = C150.subs(W321NF_eval)
    
    # C114 = C114.subs(W13NF_eval)
    # C114 = C114.subs(W33NF_eval)
    # C114 = C114.subs(W121NF_eval)
    # C114 = C114.subs(W1221NF_eval)
    # C114 = C114.subs(W321NF_eval)
    
    # C132 = C132.subs(W13NF_eval)
    # C132 = C132.subs(W33NF_eval)
    # C132 = C132.subs(W121NF_eval)
    # C132 = C132.subs(W1221NF_eval)
    # C132 = C132.subs(W321NF_eval)
    
    # C150 = C150.subs(W02NF_eval)
    # C150 = C150.subs(W22NF_eval)
    # C150 = C150.subs(W011NF_eval)
    # C150 = C150.subs(W211NF_eval)
    
    # C114 = C114.subs(W02NF_eval)
    # C114 = C114.subs(W22NF_eval)
    # C114 = C114.subs(W011NF_eval)
    # C114 = C114.subs(W211NF_eval)
    
    # C132 = C132.subs(W02NF_eval)
    # C132 = C132.subs(W22NF_eval)
    # C132 = C132.subs(W011NF_eval)
    # C132 = C132.subs(W211NF_eval)
    
    # C150 = C150.subs(phi11NF_eval)
    # C150 = C150.subs(psi11NF_eval)
    
    # C114 = C114.subs(phi11NF_eval)
    # C114 = C114.subs(psi11NF_eval)
    
    # C132 = C132.subs(phi11NF_eval)
    # C132 = C132.subs(psi11NF_eval)
    
    file = open('Fifth-order coefficient C150.txt','w')
    file.write(latex(C150))
    file.close()
    file = open('Fifth-order coefficient C114.txt','w')
    file.write(latex(C114))
    file.close()
    file = open('Fifth-order coefficient C132.txt','w')
    file.write(latex(C132))
    file.close()
    
    print('The fifth-order coefficients were computed and saved into text files.')
    
file = open('Find codimension-two bifurcation points.txt', 'w')
file.write(str(cod2))
file.close()
    
if plot2d=='y':
    try:
        for parnum in range(2):
            parameters_on_axes[parnum] = eval(parameters_on_axes[parnum])
    except:
        print('The variable párameters_on_axes is not well defined.')
        exit()

    try:
        file=open('Initial conditions for Wave bifurcation curves.txt','w')
        
        for key in lines_to_search.keys():
            if isinstance(lines_to_search[key],list):
                for initialsolnum in range(len(lines_to_search[key])):
                    if isfloat(str(lines_to_search[key][initialsolnum])):
                        file.write(latex(eval(key)) + ',' + 
                                   latex(eval(str(lines_to_search[key][initialsolnum]))) + '\n')
            else:
                if isfloat(str(lines_to_search[key])):
                    file.write(latex(eval(key)) + ',' + latex(eval(str(lines_to_search[key]))) + '\n')
                
        file.close()    
    except:
        print('The variable lines_to_search is not well defined.')
        exit()
    
    try:
        auxpar = dict()
        
        for key in parameter_functions.keys():
            auxpar[eval(key)] = eval(parameter_functions[key])
        parameter_functions = auxpar
    except:
        parameter_functions = dict()
    
    file=open('Fixed parameter values.txt','w')
    for key in parameters.keys():
        if key not in parameters_on_axes:
            if key not in parameter_functions.keys():
                file.write(latex(key) + ',' + latex(parameters[key]) + '\n')
            else:
                file.write(latex(key) + ',' + latex(parameter_functions[key]) + '\n')
            
    file.close()
        
    file = open('Parameters on axes.txt','w')
    file.write(latex(parameters_on_axes[0]) + ',' + latex(intervalx[0]) + ',' + latex(intervalx[1]) + '\n')
    file.write(latex(parameters_on_axes[1]) + ',' + latex(intervaly[0]) + ',' + latex(intervaly[1]) + '\n')
    file.close()

    try:
        if len(names_of_parameters)==0:
            for parnum in range(2):
                names_of_parameters[parnum] = latex(parameters_on_axes[parnum])
    except:
        names_of_parameters = parameters_on_axes
        for parnum in range(2):
            names_of_parameters[parnum] = latex(names_of_parameters[parnum])
        
    file = open('Actual names of parameters.txt','w')
    file.write(names_of_parameters[0] + ',' + names_of_parameters[1])
    file.close()
    
    if not os.path.isfile('Plotter.nb'):
        shutil.copyfile(os.path.dirname(os.path.realpath(__file__)) + '\\Plotter.nb', 'Plotter.nb')
    
    print('The variables to plot the bifurcation diagram in Mathematica were correctly saved')