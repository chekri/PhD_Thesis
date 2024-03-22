This folder contains the codes used for implementng the examples in Chapter 7. To implement continuation, just run any of the given script files by running the command 'auto filname.auto'. Auto computes the solutions for the given set of parameters (set in constant files c.files) and outputs the solutions as 's.xxx', where xxx is outputfile name. Use the provided python commands to plot the bifucation plots. The folder "CTCRwithoutStability" includes four auto script files, each generates solutions with different set of paramters.

1.***Example1_surface.auto***: Continuation when $\alpha^{[2]}_{o}$ is varied for different section lengths $L_{2}$.(section 7.1.2 of thesis)

2.***Example1_surface_load.auto***: Continuation when $\alpha^{[2]}_{o}$ is varied for different section lengths $L_{2}$ with a tip load.(section 7.1.3 of thesis) 

3.***Example2_surface.auto***: Continuation when $\alpha^{[2]}_{o}$ and  $\alpha^{[3]}_{o}$ is varied for different section lengths in a three-tube CTCR.(section 7.2.1 of thesis)

4.***Example3_surface.auto***: Continuation when rotating tip load acts on a CTCR.(section 7.2.2 of thesis)


The folder "CTCRwithJacobi_equations" includes 4 script files

1.***Example1.auto***: (section 7.1.2 of thesis)

2.***Example1_load.auto***: (7.1.3 of thesis) 

3.***Example2.auto***: (section 7.2.1 of thesis)

4.***Example3.auto***: (section 7.2.2 of thesis)

These scripts generate solutions conssiting of both CTCR equilibria and the corresponding Jacobi equations. The Python code "CTCR_conjugatepoints.py" computes the conjugate points from this solution files and plots the conjugate point computations. 





