function [Jac,Jac1]=Jac_Sweep_area_test(Tr_Current,Tr_Prev,U_Current,U_Prev)
PT=Sweep_area_test(Tr_Current,Tr_Prev,U_Current,U_Prev)
e=1e-02;
diff=[[e,0,0,0,0,0];[0,e,0,0,0,0];[0,0,e,0,0,0];[0,0,0,e,0,0];[0,0,0,0,e,0];[0,0,0,0,0,e]];
for i=[1:6]
    UU=U_Current+diff(i,:);
    Tr_UU=Trajectory(UU,[0,0,0]);
    Tp_f=Sweep_area_test(Tr_UU,Tr_Prev,UU,U_Prev);
    Jac(i)=(Tp_f - PT)/e;  
end

for i=[1:6]
    UPv=U_Prev+diff(i,:);
    Tr_UPv=Trajectory(UPv,[0,0,0]);
    Tp_f=Sweep_area_test(Tr_Current,Tr_UPv,U_Current,UPv);
    Jac1(i)=(Tp_f - PT)/e;  
end