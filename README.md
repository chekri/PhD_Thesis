# Instructions for Using Codes and Tools in my PhD Thesis:
These are the codes used in the examples of my PhD thesis. These programs are written in AUTO-07p, Python3 and Matlab-R2023b. The examples in part 1 and Part 2 of the thesis pertaining to numerical continuation are performed using AUTO-07p software package. 

## AUTO-07p Configuration

AUTO-07p is an easy to use software package. To install this software package, unzip the provided "auto-07p-master.zip file" and navigate to its parent folder. Then configure the files by running '.\configure' in the terminal. After successful configuration, execute 'make all'. Then, add the path of the bin folder to your system path file file by appending the line "path_to_auto/bin" to your bash file. The installation can be verified by running 'auto' in the terminal. Once installed, the codes can be run easily through auto scripts (auto.filename). To run these scripts, run auto auto.filename in the terminal window. These scripts make use of input files (`filename.f90`) and constants files (`c.constants`) where the parameters running the continuation can be tuned. The variables in c.constants files and script files can be adjusted as per the requirement. Once, the computation is complete, the solution files are written in a "s.filename" file. An example of soluton file (s.file) can be found in "Kirchoff_Rod/Conserved_quantites".  These files are analyzed and plotted using python tools provided in the respective folders. Ensure that the necessary Python dependencies, including numpy and matplotlib are installed.

## Instructions for using the codes in the folder "Kirchoff_Rod": 

### Bifurcation Surface Generation

There are three types of codes in this folder. First type of code is used for generating bifurcation surfaces. Navigate to the folder "Kirchhoff_Rod" and execute the desired script file. Six script files corresponding to six cases can be seen here, namely 1."Rod_FixedL_varyingF.auto", 2."Rod_FixedF_varyingL.auto", 3."Rod_arm_fixedL_varyingF.auto", 4."Rod_arm_fixedF_varyingL.auto", 5.Rod_torsion_fixedL_varyingF.auto, 6."Rod_torsion_fixedF_varyingL.auto".

Run the  script file, say Rod_FixedL_varyingF.auto by typing "auto Rod_FixedL_varyingF.auto" in terminal. As the name suggests, this script file generates numerical continuation solutions in $\theta$ for different Tip Loads $F$ at fixed Length $L$. Use the Python plotting tools in the "Plotting_tools" folder to generate bifurcation surface plots. Run the python code followed by the solution file (s.filename) to obtain the corresponding bifurcation surface plots. Make sure that you are running the python plotting code with appropriate solution files. The python code "surface_plot_fixed_length.py" with tag "fixed_length" works only with the solution files generated with the cases of fixed length such as 1,3,5. Otherwise, the program terminates with an error. Therefore, use solutions generated from 1,3,5 with "surface_plot_fixed_length.py" and solutions from 2,4,6
with "surface_plot_fixed_load.py" to generate the respective bifurcation surfaces.

### Stability Analysis

Second type of code is used for analyzing the stability of equilibrium by computing conjugate points. To perform these computations, access the folder "Kirchhoff_rod_with_Jacobi_equations" and run "auto data_for_conjugatepoints.auto" in the terminal. This run will generate the solutions pertaining to equilibria along with the solutions of Jacobi equations in the solution file (s.filename). Then, run the python code "Conjugate_points.py" folowed by the solution file (s.filename) to visualize the computation of conjugate points. The labels of the equilibria for which the conjugate point computations are displayed can be controlled inside the python file. These labels can be acccessed during the auto run or also can be set in the constants file (c.filename) of auto. Different cases of equilibria can be analyzed by controlling the parameters in the script file "auto data_for_conjugatepoints.auto". Jupyter notebook version for these comutations is also provided.

### Benchmarks for Numerical methods 

Third type of code is used for verifying the numerical behavior of the solutions by computing the conserved quantites of the Kirchhoff Rod system. Proceed to the folder named "Conserved_quantites" and run "auto data_for_conserved_quantites.auto" from terminal. It will generate a continuation solution file corresponding to the numerical parameters specified in the script file. Run the python code Integrals.py followed by the generated data file to plot the conserved quantites along the solutions. The numerical parameters such as mesh size "NTST" and number of collocation points "NCOL" can be varied in the script file. Different cases of solutions can also be considered by varying the parameters inside it.\\

## Instructions for using the codes in the folder "CTCR":

### Bifurcation Surface Generation

The implementation process for CTCR codes mirrors that of Kirchhoff Rod codes, with similar steps for bifurcation surface generation, stability analysis, and numerical behavior verification. Here, also we have three kinds of codes. First kind of codes is used for generating the bifurcation surfaces depicted in the thesis. Proceed to the folder "CTCR/CTCR_withoutStability" and run the script file by running the command "auto auto.filename". There are three script files in this folder corresponding to the three examples in the Part 2 of the thesis. After generating the solution files, surfaces plots can generated by the python plotting tools in the folder "CTCR/CTCR_withoutStability/Plotting_Tools". There are three plotting tools in this folder. 

There are three Python plotting tools available, each designed to work with specific solution files:

1. **Example1_surface_plot.py**: Use this script with solution files generated by scripts tagged as "Example1" (e.g., Example1_surface.auto).
2. **Example2_surface_plot.py**: Compatible with solution files from scripts labeled as "Example2" (e.g., Example2_surface.auto).
3. **Example3_surface_plot.py**: Intended for use with solution files generated by scripts tagged as "Example3" (e.g., Example3_surface.auto).

Ensure you run the appropriate Python plotting tool corresponding to the solution files you have generated. Using an incorrect tool with mismatched solution files may result in errors. Different scenarios of loadings can be considered by varying the parameters in auto script files or constants file (c.filename).

### Stability Analysis

Second type of code pertains to the stability analysis of the CTCR equilibria and the implementation is similar to that of the "Kirchhoff_Rod". Access the folder "CTCR/CTCRwith_Jacobiequations" and run any of the given script files by running "auto auto.filename" in the terminal. Replace "filename" with the actual name of the file name in the folder. This run will generate the solutions along with the solutoins of the Jacobi equations in the solution file  s.<filename>. Then, execute the Python code by runnijng "python3 CTCR_conjugate_points.py s.<filename>". This will display the conjugate point computations along the CTCR sections for the given labelled equilibrium from the solution file. This label can be varied in inside the code "CTCR_conjugate_points.py" (almost at the end of the code). Different cases of running paramters can be considered by varying parameters in scripts files (auto.filename) and constants file (c.constants).

### Benchmarks for Numerical methods 

Third type of code is used for verifying the numerical behavior of the solutions by computing the conserved quantites of the CTCR system. Proceed to the folder named "CTCR/CTCRwithut_Stability/Conserved_qunatites" and run "auto data_for_conserved_quantites.auto" from terminal. It will generate a solution file corresponding to the continuation parameters specified in the script file. Run the python code "Conserved.py" followed by the generated data file to plot the conserved quantites along the solutions. The plots along the three sections are displayed in different subplots. The numerical parameters such as mesh size "NTST" and number of collocation points "NCOL" can be varied in the script file. Different cases of solutions can also be considered by varying the parameters inside it.

## Instructions for using the codes in the folder "Optimal Control of CTCR":

The optimiation codes are implemented in MatlabR2023b.This folder contains codes correspoding to three kind of examples discussed in the chapter 10 of the dissertation, namely

1. **Example1**: Examples pertaining to the section "Minimimum Working volume"
2. **Example2**: Examples pertaining to the section "Fixed Tip Orientation with adjustable Position".
3. **Example3**: Examples pertaining to the section "Fixed Tip Position with adjustable Orientation".

Naviagte to any folder and open the .m file which starts with "optimize", such as "Optimize_path_alongpath.m" or "Optimize_path_sweep_area.m" , Then, run the command with with Optimize_path_alongpath([]) to start the optimization implementation. It can also be started from a solution from previous run. In that case, run "Optimize_path_alongpath(prev_solution)". The parameters, penalization weights, tolerances can be set and varied inside this code. Once, the computation is successful, run "Animate(output)" to postprocess the results. The evolution of CTCR configurations during the optimization process, evolution of control paramters and some relevant functions are displyed in this run.  
