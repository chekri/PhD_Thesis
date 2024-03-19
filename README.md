$\textbf{PhD Thesis}$
These are the codes used in the examples of my PhD thesis. These programs are written in AUTO-07p, python and Matlab. The examples in part 1 and Part 2 of the thesis pertaining to numerical continuation are performed using AUTO-07p software package. AUTO-07p is an easy to use software package. To install this software package, unzip the provided auto-07p-master.zip file and move to its parent folder. Then configure the files by .\configure . After the successful configuration, 'make all'. Once installed, the codes can be run easily through  auto scripts (auto.filename). To run these scripts, run auto auto.filename in the terminal window. Once, the computation is complete, the solution files are written in a "s.filename" file. These files are analyzed and plotted using python tools provided in the respective folders.

$\textbf{Instructions for using the codes in the folder} \mathbf{"Kirchoff_Rod"$: 
There are three types of codes in this folder. First type of code is used for generating bifurcation surfaces. Go inside the folder named "Kirchhoff_Rod" and run a script file. Six script files corresponding to six cases can be seen here, namely 1."Rod_FixedL_varyingF.auto", 2."Rod_FixedF_varyingL.auto", 3."Rod_arm_fixedL_varyingF.auto", 4."Rod_arm_fixedF_varyingL.auto", 5.Rod_torsion_fixedL_varyingF.auto, 6."Rod_torsion_fixedF_varyingL.auto".

Run the  script file, say Rod_FixedL_varyingF.auto by typing "auto Rod_FixedL_varyingF.auto" in terminal. As the name suggests, this script file generates numerical continuation solutions in $\theta$ for different Tip Loads $F$ at fixed Length $L$. Bifurcation surfaces can be generated by running the python codes in the folder named "Plotting_tools". Run the python code followed by the solution file (s.filename) to obtain the corresponding bifurcation surface plots. Make sure that you are running the python plotting code with appropriate solution files. The python code "surface_plot_fixed_length.py" with tag "fixed_length" works only with the solution files generated with the cases of fixed length such as 1,3,5. Otherwise, the program terminates with an error. Therefore, use solutions generated from 1,3,5 with "surface_plot_fixed_length.py" and solutions from 2,4,6
with "surface_plot_fixed_load.py" to generate the respective bifurcation surfaces.
 
Second type of code is used for analyzing the stability of equilibrium by computing conjugate points. To perform these computations, go to the folder, "Kirchhoff_rod_with_Jacobi_equations" and run "auto data_for_conjugatepoints.auto" in the terminal. This run will generate the solutions pertaining to equilibria along with the solutions of Jacobi equations. Then, run the python code "Conjugate_points.py" folowed by the generated solution file to visualize the conjugate points. The labels of the equilibria for which the conjugate point computations are displayed can be controlled inside the python file. Different cases of equilibria can be analyzed by controlling the parameters in the script file "auto data_for_conjugatepoints.auto". Jupyter notebook version for these comutations is also provided.

Third type of code is used for verifying the numerical behavior of the solutions by computing the conserved quantites of the Kirchhoff Rod system. Go to the folder named "Conserved_quantites" and run "auto data_for_conserved_quantites.auto" from terminal. It will generate a solution file corresponding to the continuation parameters specified in the script file. Run the python code Integrals.py followed by the generated data file to plot the conserved quantites along the solutions. The numerical parameters such as mesh size "NTST" and number of collocation points "NCOL" can be varied in the script file. Different cases of solutions can also be considered by varying the parameters inside it.

$\textbf{Instructions for using the codes in the folder "CTCR"}$:
The implementation of these codes is similar to that of the codes in "Kirchhoff_Rod". Here, also we have three kinds of codes. First kind of codes is used for generating the bifurcation surfaces depicted in the thesis. Go to the folder "CTCR/CTCR_withoutStability" and run the script file by running the command "auto auto.<filename>". There are three script files in this folder corresponding to the three examples in the Part 2 of the thesis. After generating the solution files, surfaces plots can generated by the python plotting tools in the folder "CTCR/CTCR_withoutStability/Plotting_Tools". There are three plotting tools in this folder. These python file with tag "Example1" works only with the solution files generated by script files with tag "Example1" and so on. For example, the python "Example1_surface_plot.py" code generates plots only with the solution files (s.<filename>) generated by the script "Example1_surface.auto" or "Example1_surface_load.auto" and generates errors with other solution files. Same applies with other remaining two python codes. Different scenarios of loadings can be considered by varying the parameters in auto script files or constants file (c.filename).

Second type of code pertains to the stability analysis of the CTCR equilibria and the implementation is similar to that of the "Kirchhoff_Rod". Go the folder "CTCR/CTCRwith_Jacobiequations" and run any of the given script files by running "auto <auto.scriptname>" in the terminal. This run will generate the solutions along with the solutoins of the Jacobi equations in the solution file  s.<filename>. Then, run the python plotting tool by running "python3 CTCR_conjugate_points.py s.<filename>". This will display the conjugate point computations along the CTCR sections for the given labelled equilibrium from the solution file. This label can be varied in inside the code "CTCR_conjugate_points.py" (almost at the end of the code). Different cases of running paramters can be considered by varying parameters in scripts files (auto.<scriptname>) and constants file (c.name).

Third type of code is used for verifying the numerical behavior of the solutions by computing the conserved quantites of the CTCR system. Go to the folder named "CTCR/CTCRwithut_Stability/Conserved_qunatites" and run "auto data_for_conserved_quantites.auto" from terminal. It will generate a solution file corresponding to the continuation parameters specified in the script file. Run the python code "Conserved.py" followed by the generated data file to plot the conserved quantites along the solutions. The plots along the three sections are displayed in different subplots. The numerical parameters such as mesh size "NTST" and number of collocation points "NCOL" can be varied in the script file. Different cases of solutions can also be considered by varying the parameters inside it.
