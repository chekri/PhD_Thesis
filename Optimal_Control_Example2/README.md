This folder contains the codes used for optimization in Example 2 (Section: "Fixed Tip Orientation with Adjustable Position") in Chapter 9 of the thesis. The CTCR maneuvers such that its tip moves to the target point along a prescribed path (straight line) and its orientation is constrained. The orientation is constrained either by penalizing its deviation or using path constraints.

1.***Example2a_inverted***: CTCR maneuver with its orientation constrained by penalizing the deviation.

2.***Example2b_inverted***: CTCR maneuver with its orientation constrained by path constraints.

3.***Example2c***: CTCR maneuvers to a farther target point with its orientation constrained by path constraints.

Navigate to any folder and open the "Optimize_follow_path.m" file and run "Optimize_follow_path([])" in the command window to start the Optimization task from an initial solution, with stationary solution (i.e., initial control parameters at all time steps). The run can also be started from the known solution "initial_solution" by running "Optimize_follow_path(initial_solution)", and in this case, the solution provided should have the same mesh size. Once the run is complete, implement "Animate(output)" to post-process the results from "output".
