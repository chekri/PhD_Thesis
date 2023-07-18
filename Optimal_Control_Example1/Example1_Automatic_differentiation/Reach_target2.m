function [Rt,Jac]=Reach_target2(Utar,U)
%This "Target-pursuit" quantifies the final the distance between the end
%point and the specified target point.
%Target point


%Distance between Target point and the current end point.
Rt=norm(Utar-U);
Jac=(U-Utar)/Rt;