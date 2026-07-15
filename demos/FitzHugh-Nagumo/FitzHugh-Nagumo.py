parameters={"a2": 0.4497,
        "a3": 0.5,
        "a0": -0.1,
        "epsilon2": 0.2,
        "epsilon3": 1,
        "D1": 0.005,
        "D2": 0,
        "D3": 1}

var=['u',
     'v',
     'w']

diffmatrix=[['D1', 0, 0],
            [0, 'D2', 0],
            [0, 0, 'D3']]

kinetics=["u - u**3 - v",
        "epsilon2*(u - a2*v - a3*w - a0)",
        "epsilon3*(u - w)"]

tol=1e-7

phiunit='n'

crosscoef='y'

crosspar='a2'

equilibrium=[]

fifthcoef='y'

orthogonal='n'

cod2='y'

plot2d='y'

parameters_on_axes=['a2', 'D2']

names_of_parameters=['a_v', 'D_2']

intervalx = [0.42, 0.436]

intervaly = [0, 0.3]

lines_to_search={'D2':0.2}

plot3d='n'

parameters_on_axes3=[]

name_of_extra_parameter=''

extrainterval=[]

extraturing=[]