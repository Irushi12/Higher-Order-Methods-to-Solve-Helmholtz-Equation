%% Generates the table of testing order of convergence with respect to maximum norm

clear all;
u1=1; % boundary condition at x=0
uend=-1; %boundary condition at x=pi
f=@(x)  3599*cos(x);% function on RHS 3599*cos(x);
k = 3600; %function of k^2(x) 3600
N=11; % number of intervals+1, so i=1,2,...N

ue1=u_exact(N)';

u_one1 = cm_4_dm(u1,uend,f,k,N);
u_two1 = cm_4_cg(u1,uend,f,k,N);
u_three1 = ccm_4_dm(u1,uend,f,k,N);
u_four1 = ccm_4_cg(u1,uend,f,k,N);
u_five1 = cm_6_dm(u1,uend,f,k,N);
u_six1 = cm_6_cg(u1,uend,f,k,N);
u_seven1 = ccm_6_dm(u1,uend,f,k,N);
u_eight1 = ccm_6_cg(u1,uend,f,k,N);
u_nine1 = cm_6_bc4_dm(u1,uend,f,k,N);
u_ten1 = cm_6_bc4_cg(u1,uend,f,k,N);
u_eleven1 = ccm_6_bc4_dm(u1,uend,f,k,N);
u_twelve1 = ccm_6_bc4_cg(u1,uend,f,k,N);

N=21;
ue2=u_exact(N)';

u_one2 = cm_4_dm(u1,uend,f,k,N);
u_two2 = cm_4_cg(u1,uend,f,k,N);
u_three2 = ccm_4_dm(u1,uend,f,k,N);
u_four2 = ccm_4_cg(u1,uend,f,k,N);
u_five2 = cm_6_dm(u1,uend,f,k,N);
u_six2 = cm_6_cg(u1,uend,f,k,N);
u_seven2 = ccm_6_dm(u1,uend,f,k,N);
u_eight2 = ccm_6_cg(u1,uend,f,k,N);
u_nine2 = cm_6_bc4_dm(u1,uend,f,k,N);
u_ten2 = cm_6_bc4_cg(u1,uend,f,k,N);
u_eleven2 = ccm_6_bc4_dm(u1,uend,f,k,N);
u_twelve2 = ccm_6_bc4_cg(u1,uend,f,k,N);

N=41;
ue3=u_exact(N)';

u_one3 = cm_4_dm(u1,uend,f,k,N);
u_two3 = cm_4_cg(u1,uend,f,k,N);
u_three3 = ccm_4_dm(u1,uend,f,k,N);
u_four3 = ccm_4_cg(u1,uend,f,k,N);
u_five3 = cm_6_dm(u1,uend,f,k,N);
u_six3 = cm_6_cg(u1,uend,f,k,N);
u_seven3 = ccm_6_dm(u1,uend,f,k,N);
u_eight3 = ccm_6_cg(u1,uend,f,k,N);
u_nine3 = cm_6_bc4_dm(u1,uend,f,k,N);
u_ten3 = cm_6_bc4_cg(u1,uend,f,k,N);
u_eleven3 = ccm_6_bc4_dm(u1,uend,f,k,N);
u_twelve3 = ccm_6_bc4_cg(u1,uend,f,k,N);

x1=log2(max(abs(u_one1-ue1))/max(abs(u_one2-ue2)));
x2=log2(max(abs(u_two1-ue1))/max(abs(u_two2-ue2)));
x3=log2(max(abs(u_three1-ue1))/max(abs(u_three2-ue2)));
x4=log2(max(abs(u_four1-ue1))/max(abs(u_four2-ue2)));
x5=log2(max(abs(u_five1-ue1))/max(abs(u_five2-ue2)));
x6=log2(max(abs(u_six1-ue1))/max(abs(u_six2-ue2)));
x7=log2(max(abs(u_seven1-ue1))/max(abs(u_seven2-ue2)));
x8=log2(max(abs(u_eight1-ue1))/max(abs(u_eight2-ue2)));
x9=log2(max(abs(u_nine1-ue1))/max(abs(u_nine2-ue2)));
x10=log2(max(abs(u_ten1-ue1))/max(abs(u_ten2-ue2)));
x11=log2(max(abs(u_eleven1-ue1))/max(abs(u_eleven2-ue2)));
x12=log2(max(abs(u_twelve1-ue1))/max(abs(u_twelve2-ue2)));

y1=log2(max(abs(u_one2-ue2))/max(abs(u_one3-ue3)));
y2=log2(max(abs(u_two2-ue2))/max(abs(u_two3-ue3)));
y3=log2(max(abs(u_three2-ue2))/max(abs(u_three3-ue3)));
y4=log2(max(abs(u_four2-ue2))/max(abs(u_four3-ue3)));
y5=log2(max(abs(u_five2-ue2))/max(abs(u_five3-ue3)));
y6=log2(max(abs(u_six2-ue2))/max(abs(u_six3-ue3)));
y7=log2(max(abs(u_seven2-ue2))/max(abs(u_seven3-ue3)));
y8=log2(max(abs(u_eight2-ue2))/max(abs(u_eight3-ue3)));
y9=log2(max(abs(u_nine2-ue2))/max(abs(u_nine3-ue3)));
y10=log2(max(abs(u_ten2-ue2))/max(abs(u_ten3-ue3)));
y11=log2(max(abs(u_eleven2-ue2))/max(abs(u_eleven3-ue3)));
y12=log2(max(abs(u_twelve2-ue2))/max(abs(u_twelve3-ue3)));
disp("------------------------------------------------------------------------")
disp("Method             h       maximum_error            order_of_convergence")
disp("------------------------------------------------------------------------")
fprintf("CM-4-DM          \pi/10         %e \n",max(abs(u_one1-ue1)))
fprintf("                 \pi/20         %e \t\t%f \n",max(abs(u_one2-ue2)),x1)
fprintf("                 \pi/40         %e \t\t%f \n",max(abs(u_one3-ue3)),y1)
fprintf("CCM-4-DM         \pi/10         %e \n",max(abs(u_three1-ue1)))
fprintf("                 \pi/20         %e \t\t%f \n",max(abs(u_three2-ue2)),x3)
fprintf("                 \pi/40         %e \t\t%f \n",max(abs(u_three3-ue3)),y3)
fprintf("CM-6-DM          \pi/10         %e \n",max(abs(u_five1-ue1)))
fprintf("                 \pi/20         %e \t\t%f \n",max(abs(u_five2-ue2)),x5)
fprintf("                 \pi/40         %e \t\t%f \n",max(abs(u_five3-ue3)),y5)
fprintf("CCM-6-DM         \pi/10         %e \n",max(abs(u_seven1-ue1)))
fprintf("                 \pi/20         %e \t\t%f \n",max(abs(u_seven2-ue2)),x7)
fprintf("                 \pi/40         %e \t\t%f \n",max(abs(u_seven3-ue3)),y7)
fprintf("CM-6-BC4-DM      \pi/10         %e \n",max(abs(u_nine1-ue1)))
fprintf("                 \pi/20         %e \t\t%f \n",max(abs(u_nine2-ue2)),x9)
fprintf("                 \pi/40         %e \t\t%f \n",max(abs(u_nine3-ue3)),y9)
fprintf("CCM-6-BC4-DM     \pi/10         %e \n",max(abs(u_eleven1-ue1)))
fprintf("                 \pi/20         %e \t\t%f \n",max(abs(u_eleven2-ue2)),x11)
fprintf("                 \pi/40         %e \t\t%f \n",max(abs(u_eleven3-ue3)),y11)
fprintf("CM-4-CG          \pi/10         %e \n",max(abs(u_two1-ue1)))
fprintf("                 \pi/20         %e \t\t%f \n",max(abs(u_two2-ue2)),x2)
fprintf("                 \pi/40         %e \t\t%f \n",max(abs(u_two3-ue3)),y2)
fprintf("CCM-4-CG         \pi/10         %e \n",max(abs(u_four1-ue1)))
fprintf("                 \pi/20         %e \t\t%f \n",max(abs(u_four2-ue2)),x4)
fprintf("                 \pi/40         %e \t\t%f \n",max(abs(u_four3-ue3)),y4)
fprintf("CM-6-CG          \pi/10         %e \n",max(abs(u_six1-ue1)))
fprintf("                 \pi/20         %e \t\t%f \n",max(abs(u_six2-ue2)),x6)
fprintf("                 \pi/40         %e \t\t%f \n",max(abs(u_six3-ue3)),y6)
fprintf("CCM-6-CG         \pi/10         %e \n",max(abs(u_eight1-ue1)))
fprintf("                 \pi/20         %e \t\t%f \n",max(abs(u_eight2-ue2)),x8)
fprintf("                 \pi/40         %e \t\t%f \n",max(abs(u_eight3-ue3)),y8)
fprintf("CM-6-BC4-CG      \pi/10         %e \n",max(abs(u_ten1-ue1)))
fprintf("                 \pi/20         %e \t\t%f \n",max(abs(u_ten2-ue2)),x10)
fprintf("                 \pi/40         %e \t\t%f \n",max(abs(u_ten3-ue3)),y10)
fprintf("CCM-6-BC4-CG     \pi/10         %e \n",max(abs(u_twelve1-ue1)))
fprintf("                 \pi/20         %e \t\t%f \n",max(abs(u_twelve2-ue2)),x12)
fprintf("                 \pi/40         %e \t\t%f \n",max(abs(u_twelve3-ue3)),y12)
