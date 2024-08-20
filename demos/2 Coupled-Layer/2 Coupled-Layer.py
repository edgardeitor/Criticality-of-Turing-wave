parameters={"epsilon1": 0.215,
        "f1": 1.1,
        "q1": 0.01,
        "delta1": 0.43,
        "epsilon2": 0.5,
        "f2": 0.65,
        "q2": 0.01,
        "delta2": 1.0,
        "D1": 0.1,
        "D2": 0.1,
        "D3": 0.1,
        "D4": 3.0,
        "D5": 100.0}

var=['u',
     'v',
     'w',
     'y',
     'z']

diffmatrix=[['D1', 0, 0, 0, 0],
            [0, 'D2', 0, 0, 0],
            [0, 0, 'D3', 0, 0],
            [0, 0, 0, 'D4', 0],
            [0, 0, 0, 0, 'D5']]

kinetics=["(u - u**2 - f1*v*(u - q1)/(u + q1))/epsilon1 - (u - w)/delta1",
        "u - v",
        "(u - w)/delta1 + (y - w)/delta2",
        "(y - y**2 - f2*z*(y - q2)/(y + q2))/epsilon2 - (y - w)/delta2",
        "y - z"]

tol=1e-7

parameter_functions={}

phiunit='n'

crosscoef='n'

crosspar=''

equilibrium=[]

fifthcoef='y'

orthogonal='n'

cod2='y'

plot2d='y'

parameters_on_axes=['f1', 'D4']

names_of_parameters = []

intervalx = [0.5, 2.5]

intervaly = [0.0, 10.0]

lines_to_search = {'f1': [1, 2]}

plot3d='n'

parameters_on_axes3=[]

name_of_extra_parameter=''

extrainterval=[]

extraturing=[]