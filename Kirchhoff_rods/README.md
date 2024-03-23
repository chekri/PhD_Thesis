This folder contains the codes used for implementing the examples in Chapter 5. To implement continuation, just run any of the given script files by running the command 'auto filename.auto'. Auto computes the solutions for the given set of parameters (set in constant files c.files) and outputs the solutions as 's.xxx', where xxx is the output file name. Use the provided Python commands to plot the bifurcation plots. The folder "Kirchoff_rods" includes six auto script files, each generating solutions with a different set of parameters.

1.***Rod_FixedL_varyingF.auto***: The bifurcation results when the intrinsically curved rod is rotated with zero arm, zero torsion, fixed Length, and varying tip load. (section 5.3.1)

2.***Rod_FixedF_varyingL.auto***: The bifurcation results when the intrinsically curved rod is rotated with zero arm, zero torsion, fixed tip Load and varying length. (section 5.3.2)

3.***Rod_arm_FixedL_varyingF.auto***: The bifurcation results when the intrinsically curved rod is rotated with a non-zero arm, zero torsion, fixed Length and varying tip load. (section 5.3.3)

4.***Rod_arm_FixedF_varyingL.auto***: The bifurcation results when the intrinsically curved rod is rotated with a non-zero arm, zero torsion, fixed tip Load and varying length. (section 5.3.3)

5.***Rod_torsion_FixedL_varyingF.auto***: The bifurcation results when the intrinsically curved rod is rotated with zero arm, non-zero torsion, fixed length and varying tip load. (section 5.3.3)

6.***Rod_torsion_FixedF_varyingL.auto***: The bifurcation results when the intrinsically curved rod is rotated with zero arm, non-zero torsion, fixed tip load and varying length. (section 5.3.3)


The folder "Kircchoff_rods/Kirchhoff_rod_with_Jacobi_equations" includes a script file, which generates both rod equilibria and the corresponding Jacobi equations. The Python code "Conjugate_points.py" computes the conjugate points from this solution files and plots the conjugate point computations. 

The folder "Kirchhoff_rods/Conserved_quantities" contains codes for verifying the convergence of the numerical solutions (section 5.4). The auto script files generate the solutions. The Python code "Conserved.py" computes the conserved quantities such as Hamiltonian, the norm of the quaternion, and $\mu \cdot \mathbf{q} + 2 \mathbf{r} \cdot \mathbf{n} $ along the length of the rod.



