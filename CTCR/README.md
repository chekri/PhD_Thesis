This folder contains the codes used for implementing the examples in Chapter 7. To implement continuation, just run any of the given script files by running the command 'auto filename.auto'. Auto computes the solutions for the given set of parameters (set in constant files c.files) and outputs the solutions as 's.xxx', where xxx is the output file name. Use the provided Python commands to plot the bifurcation plots. The folder "CTCRwithoutStability" includes four auto script files, each generating solutions with a different set of parameters.

1.***Example1_surface.auto***: Continuation when $\alpha^{[2]}_{o}$ is varied for different section lengths $L_{2}$.(section 7.1.2 of thesis)

2.***Example1_surface_load.auto***: Continuation when $\alpha^{[2]}_{o}$ is varied for different section lengths $L_{2}$ with a tip load.(section 7.1.3 of thesis) 

3.***Example2_surface.auto***: Continuation when $\alpha^{[2]}_{o}$ and  $\alpha^{[3]}_{o}$ is varied for different section lengths in a three-tube CTCR.(section 7.2.1 of thesis)

4.***Example3_surface.auto***: Continuation when rotating tip load acts on a CTCR.(section 7.2.2 of thesis)


The folder "CTCRwithJacobi_equations" includes 4 script files

1.***Example1.auto***: (section 7.1.2 of thesis)

2.***Example1_load.auto***: (section 7.1.3 of thesis) 

3.***Example2.auto***: (section 7.2.1 of thesis)

4.***Example3.auto***: (section 7.2.2 of thesis)

These scripts generate solutions consisting of both CTCR equilibria and the corresponding Jacobi equations. The Python code "CTCR_conjugatepoints.py" computes the conjugate points from this solution files and plots the conjugate point computations. 

The folder "CTCRwithoutStability/Conserved_quantities" contains codes for verifying the convergence of the numerical solutions (section 7.3). The auto script files generate the solutions. The Python code "Conserved.py" computes the conserved quantities such as Hamiltonian, the norm of the quaternion, and $\mu \cdot \mathbf{q} + 2 \mathbf{r} \cdot \mathbf{n} $ along different sections of the tubes.





