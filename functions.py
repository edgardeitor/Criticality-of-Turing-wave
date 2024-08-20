def isfloat(value):
  try:
    float(eval(value))
    return True
  except:
    return False

class Vector:
    def __init__(self,name):
        self.dummy = []
        self.actualcoord = []
        for varnum in range(nvar):
            self.dummy.append(symbols(name + '_' + str(varnum)))
            self.actualcoord.append(symbols(name + '_' + str(varnum)))
        self.dummy = Matrix(self.dummy)
        self.conj = self.conjugatevector(self.dummy)
    class conjugatevector:
        def __init__(self,vector):
            self.dummy = conjugate(vector)
            
class matrix:
    def __init__(self, name, ncol = nvar, refmatrix = ones(nvar)):
        self.dummy = zeros(ncol)
        if refmatrix!=ones(nvar):
            self.actualcoord = refmatrix
        for row in range(ncol):
            for col in range(ncol):
                if refmatrix[row,col]!=0:
                    self.dummy[row,col] = symbols(name + '_' + str(row+1) + '^' + str(col+1))
                    
def dummyvareval(vector, negativeRHS, coefmat):
    for row in range(nvar):
        vector = vector.subs(negativeRHS.dummy[row], negativeRHS.actualcoord[row])
        for col in range(nvar):
            vector = vector.subs(coefmat.dummy[row, col], coefmat.actualcoord[row, col])
    return vector

def linearsolver(vector, negativeRHS, coefmat):
    vector.actualcoord = linsolve(Add(Mul(coefmat.dummy, vector.dummy), negativeRHS.dummy), list(vector.dummy))
    vector.actualcoord = transpose(Matrix(list(vector.actualcoord)))
    vector.actualcoord = dummyvareval(vector.actualcoord, negativeRHS, coefmat)
    return vector
                    
def crossorderapplied(u, parameter, equilibrium):
    SS = zeros(nvar, 1)
    firstordereval = firstorderderivatives
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            firstordereval[counter1] = Matrix(firstordereval[counter1]).subs(var[counter2], equilibrium[counter2])
    for counter1 in range(nvar):
        SS = Add(SS, Mul(u[counter1], diff(firstordereval[counter1], parameter)))
    return SS
    
def secondorderapplied(u, v):
    DS = zeros(nvar,1)
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            DS = Add(DS, Mul(u.dummy[counter1], v.dummy[counter2],
                             secondorderderivatives[counter1][counter2]))
    return Mul(Pow(2, -1), DS)
    
def thirdorderapplied(u, v, w):
    TS = zeros(nvar,1)
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            for counter3 in range(nvar):
                TS = Add(TS, Mul(u.dummy[counter1], v.dummy[counter2], w.dummy[counter3],
                                 thirdorderderivatives[counter1][counter2][counter3]))
    return Mul(Pow(math.factorial(3), -1), TS)

def fourthorderapplied(u, v, w, r):
    Q4S = zeros(nvar,1)
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            for counter3 in range(nvar):
                for counter4 in range(nvar):
                    Q4S = Add(Q4S, Mul(u.dummy[counter1], v.dummy[counter2], w.dummy[counter3], r.dummy[counter4],
                                       fourthorderderivatives[counter1][counter2][counter3][counter4]))
    return Mul(Pow(math.factorial(4), -1), Q4S)

def fifthorderapplied(u, v, w, r, s):
    Q5S = zeros(nvar,1)
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            for counter3 in range(nvar):
                for counter4 in range(nvar):
                    for counter5 in range(nvar):
                        Q5S = Add(Q5S, Mul(u.dummy[counter1], v.dummy[counter2],
                                           w.dummy[counter3], r.dummy[counter4], s.dummy[counter5],
                                           fifthorderderivatives[counter1][counter2][counter3][counter4][counter5]))
    return Mul(Pow(math.factorial(5), -1), Q5S)

def kernel_determination(vector, coefmat, criticalcol, coefsubmatrix, submatrixrows, submatrixcols):
    vector.actualcoord[criticalcol] = 1
    
    auxiliaryterm, = linsolve(Add(Mul(coefsubmatrix.dummy, Matrix(vector.actualcoord).extract(submatrixcols, [0])),
                                  coefmat.dummy.extract(submatrixrows, [criticalcol])),
                              list(Matrix(vector.actualcoord).extract(submatrixcols, [0])))
    
    vector.actualcoord[0:criticalcol] = auxiliaryterm[0:criticalcol]
    vector.actualcoord[criticalcol+1:nvar] = auxiliaryterm[criticalcol:nvar-1]
    
    vector.actualcoord = Matrix(vector.actualcoord)
        
    for row in range(nvar):
        for col in range(nvar):
            vector.actualcoord = vector.actualcoord.subs(coefmat.dummy[row, col], coefmat.actualcoord[row, col])
            if row<nvar-1 and col<nvar-1:
                vector.actualcoord = vector.actualcoord.subs(coefsubmatrix.dummy[row, col],
                                                             coefsubmatrix.actualcoord[row, col])
    
    return vector

def critical_linearsolver(vector, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols):
    vector.actualcoord[criticalcol] = 0
    
    auxiliaryterm, = linsolve(Add(Mul(coefsubmatrix.dummy, Matrix(vector.actualcoord).extract(submatrixcols, [0])),
                                  negativeRHS.dummy.extract(submatrixrows, [0])),
                              list(Matrix(vector.actualcoord).extract(submatrixcols, [0])))
        
    vector.actualcoord[0:criticalcol] = auxiliaryterm[0:criticalcol]
    vector.actualcoord[criticalcol+1:nvar] = auxiliaryterm[criticalcol:nvar-1]
    
    vector.actualcoord = Matrix(vector.actualcoord)
    
    for row in range(nvar):
        vector.actualcoord = vector.actualcoord.subs(negativeRHS.dummy[row],negativeRHS.actualcoord[row])
        for col in range(nvar-1):
            if row<nvar-1:
                vector.actualcoord = vector.actualcoord.subs(coefsubmatrix.dummy[row,col],
                                                             coefsubmatrix.actualcoord[row,col])
                
    return vector

def evaluation_dict(vector):
    actual_dict = dict(zip(vector.dummy, vector.actualcoord))
    return actual_dict

def write_vector(vector, file):
    '''This function saves the coordinates of a vector into a text file'''
    for varnum in range(nvar):
        if varnum<nvar - 1:
            file.write(latex(vector.actualcoord[varnum]) + ',')
        else:
            file.write(latex(vector.actualcoord[varnum]) + '\n')