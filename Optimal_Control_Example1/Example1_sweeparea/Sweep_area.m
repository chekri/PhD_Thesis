function [earea,jac_t,jac_p]=Sweep_area(Tr_Current,Tr_Prev,U_Current,U_Prev,Load)
Steps=50;
e=1e-04;

IC_Prev=Tr_Prev.y(31:48,1);
IC_Current=Tr_Current.y(31:48,1);


Arc_c=[U_Current(6)*(Tr_Current.x), U_Current(6)+ U_Current(5)*(Tr_Current.x), U_Current(6)+ U_Current(5)+ U_Current(4)*(Tr_Current.x)];
Arc_p=[U_Prev(6)*(Tr_Prev.x), U_Prev(6)+ U_Prev(5)*(Tr_Prev.x), U_Prev(6)+ U_Prev(5)+ U_Prev(4)*(Tr_Prev.x)];

Cur_x=[Tr_Current.y(31,:),Tr_Current.y(1,:),Tr_Current.y(17,:)];
Cur_y=[Tr_Current.y(32,:),Tr_Current.y(2,:),Tr_Current.y(18,:)];
Cur_z=[Tr_Current.y(33,:),Tr_Current.y(3,:),Tr_Current.y(19,:)];

index=size(Tr_Current.x,2);
Arc_c(2*index)=[];
Arc_c(index+1)=[];
Cur_x(2*index)=[];
Cur_x(index+1)=[];
Cur_y(2*index)=[];
Cur_y(index+1)=[];
Cur_z(2*index)=[];
Cur_z(index+1)=[];

Prv_x=[Tr_Prev.y(31,:),Tr_Prev.y(1,:),Tr_Prev.y(17,:)];
Prv_y=[Tr_Prev.y(32,:),Tr_Prev.y(2,:),Tr_Prev.y(18,:)];
Prv_z=[Tr_Prev.y(33,:),Tr_Prev.y(3,:),Tr_Prev.y(19,:)];

index=size(Tr_Prev.x,2);
Arc_p(2*index)=[];
Arc_p(index+1)=[];

Prv_x(index+1)=[];
Prv_x(2*index)=[];
Prv_y(index+1)=[];
Prv_y(2*index)=[];
Prv_z(index+1)=[];
Prv_z(2*index)=[];

Cur_Arc=linspace(0,U_Current(4)+U_Current(5)+U_Current(6),Steps);
Prv_Arc=linspace(0,U_Prev(4)+U_Prev(5)+U_Prev(6),Steps);
%Tr_Current
%U_Current
%Arc_c
C_x=spline(Arc_c,Cur_x,Cur_Arc);
C_y=spline(Arc_c,Cur_y,Cur_Arc);
C_z=spline(Arc_c,Cur_z,Cur_Arc);

P_x=spline(Arc_p,Prv_x,Prv_Arc);
P_y=spline(Arc_p,Prv_y,Prv_Arc);
P_z=spline(Arc_p,Prv_z,Prv_Arc);
earea=enclosed_area1([C_x;C_y;C_z],[P_x;P_y;P_z]); 

E=e*[[1 0 0 0 0 0];[0 1 0 0 0 0];[0 0 1 0 0 0];[0 0 0 1 0 0];[0 0 0 0 1 0];[0 0 0 0 0 1]];
for i=[1:6]
    %U_Current+E(i,:)
    Tr_Cur=Trajectory(U_Current+E(i,:),Load,[]);
    %IC_Prev=Tr_Prev.y(31:48,1);
    %IC_Current=Tr_Current.y(31:48,1);
    Arc_c=[U_Current(6)*(Tr_Cur.x), U_Current(6)+ U_Current(5)*(Tr_Cur.x), U_Current(6)+ U_Current(5)+ U_Current(4)*(Tr_Cur.x)];
    Arc_p=[U_Prev(6)*(Tr_Prev.x), U_Prev(6)+ U_Prev(5)*(Tr_Prev.x), U_Prev(6)+ U_Prev(5)+ U_Prev(4)*(Tr_Prev.x)];

    Cur_x=[Tr_Cur.y(31,:),Tr_Cur.y(1,:),Tr_Cur.y(17,:)];
    Cur_y=[Tr_Cur.y(32,:),Tr_Cur.y(2,:),Tr_Cur.y(18,:)];
    Cur_z=[Tr_Cur.y(33,:),Tr_Cur.y(3,:),Tr_Cur.y(19,:)];

    index=size(Tr_Cur.x,2);
    Arc_c(2*index)=[];
    Arc_c(index+1)=[];
    Cur_x(2*index)=[];
    Cur_x(index+1)=[];
    Cur_y(2*index)=[];
    Cur_y(index+1)=[];
    Cur_z(2*index)=[];
    Cur_z(index+1)=[];

    Prv_x=[Tr_Prev.y(31,:),Tr_Prev.y(1,:),Tr_Prev.y(17,:)];
    Prv_y=[Tr_Prev.y(32,:),Tr_Prev.y(2,:),Tr_Prev.y(18,:)];
    Prv_z=[Tr_Prev.y(33,:),Tr_Prev.y(3,:),Tr_Prev.y(19,:)];

    index=size(Tr_Prev.x,2);
    Arc_p(2*index)=[];
    Arc_p(index+1)=[];

    Prv_x(index+1)=[];
    Prv_x(2*index)=[];
    Prv_y(index+1)=[];
    Prv_y(2*index)=[];
    Prv_z(index+1)=[];
    Prv_z(2*index)=[];

    Cur_Arc=linspace(0,U_Current(4)+U_Current(5)+U_Current(6),Steps);
    Prv_Arc=linspace(0,U_Prev(4)+U_Prev(5)+U_Prev(6),Steps);

    C_x=spline(Arc_c,Cur_x,Cur_Arc);
    C_y=spline(Arc_c,Cur_y,Cur_Arc);
    C_z=spline(Arc_c,Cur_z,Cur_Arc);

    P_x=spline(Arc_p,Prv_x,Prv_Arc);
    P_y=spline(Arc_p,Prv_y,Prv_Arc);
    P_z=spline(Arc_p,Prv_z,Prv_Arc);
    eareal=enclosed_area1([C_x;C_y;C_z],[P_x;P_y;P_z]); 
    jac_t(i)=(eareal - earea)/e;
    
end

for i=[1:6]
    Tr_Prev=Trajectory(U_Prev+E(i,:),Load,[]);
    %IC_Prev=Tr_Prev.y(31:48,1);
    %IC_Current=Tr_Current.y(31:48,1);


    Arc_c=[U_Current(6)*(Tr_Current.x), U_Current(6)+ U_Current(5)*(Tr_Current.x), U_Current(6)+ U_Current(5)+ U_Current(4)*(Tr_Current.x)];
    Arc_p=[U_Prev(6)*(Tr_Prev.x), U_Prev(6)+ U_Prev(5)*(Tr_Prev.x), U_Prev(6)+ U_Prev(5)+ U_Prev(4)*(Tr_Prev.x)];

    Cur_x=[Tr_Current.y(31,:),Tr_Current.y(1,:),Tr_Current.y(17,:)];
    Cur_y=[Tr_Current.y(32,:),Tr_Current.y(2,:),Tr_Current.y(18,:)];
    Cur_z=[Tr_Current.y(33,:),Tr_Current.y(3,:),Tr_Current.y(19,:)];

    index=size(Tr_Current.x,2);
    Arc_c(2*index)=[];
    Arc_c(index+1)=[];
    Cur_x(2*index)=[];
    Cur_x(index+1)=[];
    Cur_y(2*index)=[];
    Cur_y(index+1)=[];
    Cur_z(2*index)=[];
    Cur_z(index+1)=[];

    Prv_x=[Tr_Prev.y(31,:),Tr_Prev.y(1,:),Tr_Prev.y(17,:)];
    Prv_y=[Tr_Prev.y(32,:),Tr_Prev.y(2,:),Tr_Prev.y(18,:)];
    Prv_z=[Tr_Prev.y(33,:),Tr_Prev.y(3,:),Tr_Prev.y(19,:)];

    index=size(Tr_Prev.x,2);
    Arc_p(2*index)=[];
    Arc_p(index+1)=[];

    Prv_x(index+1)=[];
    Prv_x(2*index)=[];
    Prv_y(index+1)=[];
    Prv_y(2*index)=[];
    Prv_z(index+1)=[];
    Prv_z(2*index)=[];

    Cur_Arc=linspace(0,U_Current(4)+U_Current(5)+U_Current(6),Steps);
    Prv_Arc=linspace(0,U_Prev(4)+U_Prev(5)+U_Prev(6),Steps);

    C_x=spline(Arc_c,Cur_x,Cur_Arc);
    C_y=spline(Arc_c,Cur_y,Cur_Arc);
    C_z=spline(Arc_c,Cur_z,Cur_Arc);

    P_x=spline(Arc_p,Prv_x,Prv_Arc);
    P_y=spline(Arc_p,Prv_y,Prv_Arc);
    P_z=spline(Arc_p,Prv_z,Prv_Arc);
    eareal=enclosed_area1([C_x;C_y;C_z],[P_x;P_y;P_z]); 
    jac_p(i)=(eareal - earea)/e;
    
end

function area=enclosed_area1(R_prev,R_cur) 
area=0;
for i=[1:size(R_prev,2)-1]
    %plot3([R_cur(1,i+1),R_prev(1,i+1)],[R_cur(2,i+1),R_prev(2,i+1)],[R_cur(3,i+1),R_prev(3,i+1)],'g');
    (R_cur(:,i) - R_cur(:,i+1));
    a1=cross((R_cur(:,i) - R_cur(:,i+1)),(R_prev(:,i+1) - R_cur(:,i+1)))/2;
    a2=cross((R_prev(:,i) - R_prev(:,i+1)),(R_prev(:,i) - R_cur(:,i)))/2;
    area=area+norm(a1)+norm(a2);
end


