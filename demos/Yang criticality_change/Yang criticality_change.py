parameters={"k1": -8.5,
        "k3": 10,
        "k4": 3,
        "tau": 10,
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

parameters_on_axes=['k1','k3']

names_of_parameters=['k_0','k_v']

intervalx=[-8,-1]

intervaly=[2,10]

lines_to_search={'k3':2.5, 'k1':-2.4}

plot3d='n'

parameters_on_axes3=[]

name_of_extra_parameter=''

extrainterval=[]

extraturing=[]