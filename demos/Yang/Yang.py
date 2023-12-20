parameters={"k1": -8.5,
        "k3": 10,
        "k4": 2,
        "tau": 50,
        "D1": 1,
        "D2": 1,
        "D3": 60}

var=['u',
     'v',
     'w']

diffmatrix=[['D1', 0, 0],
            [0, 'D2', 0],
            [0, 0, 'D3']]

kinetics=["k1+2*u-u**3-k3*v-k4*w",
        "(u-v)/tau",
        "u-w"]

tol=1e-7

phiunit='n'

crosscoef='y'

crosspar='k1'

equilibrium=[]

fifthcoef='y'

orthogonal='n'

cod2='y'

plot2d='y'

parameters_on_axes=['k1','k4']

names_of_parameters=[]

intervalx=[-12,-6]

intervaly=[1,15]

lines_to_search={'k1':-8}

plot3d='n'

parameters_on_axes3=[]

name_of_extra_parameter=''

extrainterval=[]

extraturing=[]