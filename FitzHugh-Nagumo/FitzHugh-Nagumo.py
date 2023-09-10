parameters={"a_2": 0.4497,
        "a_3": 0.5,
        "a_0": -0.1,
        "epsilon_2": 0.2,
        "epsilon_3": 1,
        "D1": 0.005,
        "D2": 0,
        "D3": 1}

var=['u',
     'v',
     'w']

diffmatrix=[['D1', 0, 0],
            [0, 'D2', 0],
            [0, 0, 'D3']]

kinetics=["u-u**3-v",
        "epsilon_2*(u-a_2*v-a_3*w-a_0)",
        "epsilon_3*(u-w)"]

tol=1e-7

phiunit='n'

crosscoef='n'

crosspar=''

equilibrium=[]

fifthcoef='y'

orthogonal='n'

cod2='y'

plot2d='y'

parameters_on_axes=['a_2','D2']

names_of_parameters=['a_v','D_2']

intervalx=[0.4,0.6]

intervaly=[0,1]

lines_to_search={'D2':0.2}

plot3d='n'

parameters_on_axes3=[]

name_of_extra_parameter=''

extrainterval=[]

extraturing=[]
