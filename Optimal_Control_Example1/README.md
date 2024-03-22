This folder consists of codes used in Example 1 in the chapter 10, which model the CTCR manueveur such that its tip reaches a specified target point using minimum working volume.

1. ***Example1_Deviation_from_FTL***: Optimization uses "Deviation_from_FTL" objective and uses gradients evaluated by finite differences.
2. ***Example1_sweep_area***: Optimization uses "Sweep area" objective and uses gradients evaluated by finite differences.
 differentiation.
3. ***Example1_Automatic_differentiation***: Optimization uses "Deviation_from_FTL" objective and uses gradients evaluated by Automatic differentiation.
4. ***Example1_restricting_length_controls***: Optimization uses "Deviation_from_FTL" objective, but its length controls are restricted by using corresponding gradient components as zero.

Access any of these four folders and open the .m file that starts with "Optimize", such as Optimize_Deviation_fromFTL.m. Run the command "Optimize_Deviation_fromFTL([])" in the command window to initiate the Optimization task from an initial solution, with stationary solution (i.e., initial control parameters at all time steps). The run can also be started from the known solution "initial_solution" by running "Optimize_Deviation_fromFTL(initial_solution)", and in this case, the solution provided should have the same mesh size. Once the run is complete, implement "Animate(output)" to post-process the results from "output". 
