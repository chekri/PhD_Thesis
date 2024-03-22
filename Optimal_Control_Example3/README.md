This folder contains the codes used for optimization in Example 3 (Section: "Fixed Tip Position with Adjustable Orientation") in Chapter9 of the thesis. The CTCR maneveurs such that its tip moves to target point along a prescribeed path (straight line) and its orientation is constrained. The orientation is constrained either by penalizing its deviation or using path constraints.


1.***Example2c***: CTCR manuveur to a farther target point with its orientation constrained by path constraints.

Navigate to any folder and open the "Optimize_stay_at_point.m" file and run "Optimize_stay_at_point([])" in command window to start the Optimization task from an initial solution, with stationary solution (i.e., initial control parameters at all time steps). The run can also be started from the known solution "initial_solution" by running "Optimize_stay_at_point(initial_solution)", and in this case, the solution provided should have the same mesh size. Once the run is complete, implement "Animate(output)" to post-process the results from "output".
