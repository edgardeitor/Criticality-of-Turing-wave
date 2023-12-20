# Criticality-of-Turing-wave
This code is a general version of the code to compute the criticality of the Turing-wave bifurcation for both, travelling and standing waves.

To run the code, you must place the files 'Criticality.py', 'functions.py', and 'Plotter.nb' in the same folder, together with different folders for the models you want to study (one folder per model). The folder can be named in any way. Inside that folder, there must be a Python file 'foo.py', where 'foo' is the same name as the folder. That Python file must have the model and some standard constants to run the code. To run the code, you just need to run the file 'Criticality.py' and provide the name of the folder in which your model is stored (there will be a prompt asking for this input).

After running the Python code, there will be several new files in the folder with your model. In particular, if you get into it, you will find a file called 'Plotter.nb' that you can open with Mathematica in order to plot a bifurcation diagram coloured according to the criticality of the bifurcation and compute unfolding and fifth-order coefficients.

There are some demos in the folder "demos" that can be read and run to understand how the code works.

The constants that need to be defined for the code to run are always the same so you just need to take a demo and change its variables to analyze your models.
