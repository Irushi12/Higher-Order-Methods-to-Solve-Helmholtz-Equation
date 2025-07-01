%% Generates the table of testing order of convergence with respect to the maximun norm

clear all;
ux1 = @(x) sin(x); % boundary condition
uxend = @(x) -sin(x); %boundary condition
uy1 = @(y) 0; % boundary condition
uyend = @(y) cos(2*y); % boundary condition
k = @(x,y) 25000*(1+x^2/2+y^2/2); 
f=@(x,y) -5*sin(x)*cos(2*y)+(25000)*(1+x^2/2+y^2/2)*sin(x)*cos(2*y);% function on RHS
N=11;
ue1=u_exact2d(N);
u_1 = full(cm_4_dm_2d(N,ux1,uxend,uy1,uyend,f,k));
u_2 = full(cm_4_cg_2d(N,ux1,uxend,uy1,uyend,f,k));
u_3 = full(cm_6_dm_2d(N,ux1,uxend,uy1,uyend,f,k));
u_4 = full(cm_6_cg_2d(N,ux1,uxend,uy1,uyend,f,k));
u_5 = full(cm_6_bc4_dm_2d(N,ux1,uxend,uy1,uyend,f,k));
u_6 = full(cm_6_bc4_cg_2d(N,ux1,uxend,uy1,uyend,f,k));
u_8 = full(adi4_dm_cg(N,ux1,uxend,uy1,uyend,f,k));
u_10 = full(adi6_dm_cg(N,ux1,uxend,uy1,uyend,f,k));
u_12 = full(adi6_bc4_dm_cg(N,ux1,uxend,uy1,uyend,f,k));

N=21;
ue2=u_exact2d(N);
v_1 = full(cm_4_dm_2d(N,ux1,uxend,uy1,uyend,f,k));
v_2 = full(cm_4_cg_2d(N,ux1,uxend,uy1,uyend,f,k));
v_3 = full(cm_6_dm_2d(N,ux1,uxend,uy1,uyend,f,k));
v_4 = full(cm_6_cg_2d(N,ux1,uxend,uy1,uyend,f,k));
v_5 = full(cm_6_bc4_dm_2d(N,ux1,uxend,uy1,uyend,f,k));
v_6 = full(cm_6_bc4_cg_2d(N,ux1,uxend,uy1,uyend,f,k));
v_8 = full(adi4_dm_cg(N,ux1,uxend,uy1,uyend,f,k));
v_10 = full(adi6_dm_cg(N,ux1,uxend,uy1,uyend,f,k));
v_12 = full(adi6_bc4_dm_cg(N,ux1,uxend,uy1,uyend,f,k));


N=41;
ue3=u_exact2d(N);
w_1 = full(cm_4_dm_2d(N,ux1,uxend,uy1,uyend,f,k));
w_2 = full(cm_4_cg_2d(N,ux1,uxend,uy1,uyend,f,k));
w_3 = full(cm_6_dm_2d(N,ux1,uxend,uy1,uyend,f,k));
w_4 = full(cm_6_cg_2d(N,ux1,uxend,uy1,uyend,f,k));
w_5 = full(cm_6_bc4_dm_2d(N,ux1,uxend,uy1,uyend,f,k));
w_6 = full(cm_6_bc4_cg_2d(N,ux1,uxend,uy1,uyend,f,k));
w_8 = full(adi4_dm_cg(N,ux1,uxend,uy1,uyend,f,k));
w_10 = full(adi6_dm_cg(N,ux1,uxend,uy1,uyend,f,k));
w_12 = full(adi6_bc4_dm_cg(N,ux1,uxend,uy1,uyend,f,k));


x1=log2(max(max(abs(u_1-ue1)))/max(max(abs(v_1-ue2))));
x2=log2(max(max(abs(u_2-ue1)))/max(max(abs(v_2-ue2))));
x3=log2(max(max(abs(u_3-ue1)))/max(max(abs(v_3-ue2))));
x4=log2(max(max(abs(u_4-ue1)))/max(max(abs(v_4-ue2))));
x5=log2(max(max(abs(u_5-ue1)))/max(max(abs(v_5-ue2))));
x6=log2(max(max(abs(u_6-ue1)))/max(max(abs(v_6-ue2))));
x8=log2(max(max(abs(u_8-ue1)))/max(max(abs(v_8-ue2))));
x10=log2(max(max(abs(u_10-ue1)))/max(max(abs(v_10-ue2))));
x12=log2(max(max(abs(u_12-ue1)))/max(max(abs(v_12-ue2))));

y1=log2(max(max(abs(v_1-ue2)))/max(max(abs(w_1-ue3))));
y2=log2(max(max(abs(v_2-ue2)))/max(max(abs(w_2-ue3))));
y3=log2(max(max(abs(v_3-ue2)))/max(max(abs(w_3-ue3))));
y4=log2(max(max(abs(v_4-ue2)))/max(max(abs(w_4-ue3))));
y5=log2(max(max(abs(v_5-ue2)))/max(max(abs(w_5-ue3))));
y6=log2(max(max(abs(v_6-ue2)))/max(max(abs(w_6-ue3))));
y8=log2(max(max(abs(v_8-ue2)))/max(max(abs(w_8-ue3))));
y10=log2(max(max(abs(v_10-ue2)))/max(max(abs(w_10-ue3))));
y12=log2(max(max(abs(v_12-ue2)))/max(max(abs(w_12-ue3))));
disp("------------------------------------------------------------------------")
disp("Method             h       maximum_error            order_of_convergence")
disp("------------------------------------------------------------------------")
fprintf("CM-4-DM-2D       \pi/10         %e \n",max(max(abs(u_1-ue1))))
fprintf("                 \pi/20         %e \t\t%f \n",max(max(abs(v_1-ue2))),x1)
fprintf("                 \pi/40         %e \t\t%f \n",max(max(abs(w_1-ue3))),y1)
fprintf("CM-6-DM-2D       \pi/10         %e \n",max(max(abs(u_3-ue1))))
fprintf("                 \pi/20         %e \t\t%f \n",max(max(abs(v_3-ue2))),x3)
fprintf("                 \pi/40         %e \t\t%f \n",max(max(abs(w_3-ue3))),y3)
fprintf("CM-6-BC4-DM-2D   \pi/10         %e \n",max(max(abs(u_5-ue1))))
fprintf("                 \pi/20         %e \t\t%f \n",max(max(abs(v_5-ue2))),x5)
fprintf("                 \pi/40         %e \t\t%f \n",max(max(abs(w_5-ue3))),y5)
fprintf("CM-4-CG-2D       \pi/10         %e \n",max(max(abs(u_2-ue1))))
fprintf("                 \pi/20         %e \t\t%f \n",max(max(abs(v_2-ue2))),x2)
fprintf("                 \pi/40         %e \t\t%f \n",max(max(abs(w_2-ue3))),y2)
fprintf("CM-6-CG-2D       \pi/10         %e \n",max(max(abs(u_4-ue1))))
fprintf("                 \pi/20         %e \t\t%f \n",max(max(abs(v_4-ue2))),x4)
fprintf("                 \pi/40         %e \t\t%f \n",max(max(abs(w_4-ue3))),y4)
fprintf("CM-6-BC4-CG-2D   \pi/10         %e \n",max(max(abs(u_6-ue1))))
fprintf("                 \pi/20         %e \t\t%f \n",max(max(abs(v_6-ue2))),x6)
fprintf("                 \pi/40         %e \t\t%f \n",max(max(abs(w_6-ue3))),y6)
fprintf("ADI4-DM-CG       \pi/10         %e \n",max(max(abs(u_8-ue1))))
fprintf("                 \pi/20         %e \t\t%f \n",max(max(abs(v_8-ue2))),x2)
fprintf("                 \pi/40         %e \t\t%f \n",max(max(abs(w_8-ue3))),y2)
fprintf("ADI6-DM-CG       \pi/10         %e \n",max(max(abs(u_10-ue1))))
fprintf("                 \pi/20         %e \t\t%f \n",max(max(abs(v_10-ue2))),x4)
fprintf("                 \pi/40         %e \t\t%f \n",max(max(abs(w_10-ue3))),y4)
fprintf("ADI6-BC4-DM-CG   \pi/10         %e \n",max(max(abs(u_12-ue1))))
fprintf("                 \pi/20         %e \t\t%f \n",max(max(abs(v_12-ue2))),x6)
fprintf("                 \pi/40         %e \t\t%f \n",max(max(abs(w_12-ue3))),y6)
