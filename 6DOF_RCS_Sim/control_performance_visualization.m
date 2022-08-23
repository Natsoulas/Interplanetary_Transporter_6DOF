
%visualizations/data processing script for acs simulation
figure
plot(1:t_seconds,controllerout_log(1,:),1:t_seconds,controllerout_log(2,:),1:t_seconds,controllerout_log(3,:),'LineWidth',1)
title('Torque Command vs. Time','FontSize',12)
xlabel('Time (seconds)','FontSize',14)
ylabel('Commanded Torque (N*m)','FontSize',14)
legend('x','y','z','FontSize',12)


figure
plot(0:0.01:t_seconds-.001,Torque_res_log(1,:),0:0.01:t_seconds-.001,Torque_res_log(2,:),0:0.01:t_seconds-.001,Torque_res_log(3,:),'LineWidth',1.2)
title('Torque vs. Time','FontSize',16)
xlabel('Time (seconds)','FontSize',20)
ylabel('Torque (N*m)','FontSize',20)
legend('x','y','z','FontSize',18)

figure
plot(0:0.001:t_seconds-.001, U_s_log,'LineWidth',0.50)
title('Simplex output per thruster')
xlabel('Time (seconds)')
ylabel('Simplex output (Binary on-off)')
legend('x','y','z')


figure
plot(1:t_seconds,z_log(:,1:4),'LineWidth',0.95)
title('State Body to Inertial Quaternion vs. Time')
xlabel('Time (seconds)')
ylabel('Quaternion component value')
legend('q0','q1','q2','q3')

figure
plot(1:t_seconds,z_log(:,5:7),'LineWidth',1)
title('Angular Velocity vs. Time','FontSize',12)
xlabel('Time (seconds)','FontSize',14)
ylabel('Angular Velocity (radians/second)','FontSize',14)
legend('w1', 'w2', 'w3','FontSize',12)

figure
plot(0:0.01:t_seconds-.001,unsplinedT_log(1,:),0:0.01:t_seconds-.001,unsplinedT_log(2,:),0:0.01:t_seconds-.001,unsplinedT_log(3,:),'LineWidth',0.50)
title('Unsplined Torque vs. Time')
xlabel('Time (seconds)')
ylabel('Torque (N*m)')
legend('x','y','z')


%% 

%on and off activity horizontal lines that change color
figure
hold on 
%yyaxis left
plot(repmat([0 t_seconds],24,1)',((1:24)'.*ones(24,2))','k','LineWidth',0.90)
ylim([0,25])
yticks('manual')
yticks(0:1:25)
ylabel('Individual Thruster (number)','FontSize',14)
thing = U_s_log;thing(thing == 0) = nan;
thang = (1:24)'.*thing;
plot(repmat(0:0.001:(t_seconds-0.001),24,1)',thang','r.')
yyaxis right
rightlabelstr = strings([26,1]);
for indx = 1:size(IC.thrusterdata,1)
    rightlabelstr(indx+1) = mat2str(IC.thrusterdata(indx,5:7));
end
ylabel('Thruster firing direction','FontSize',14)
ax = gca;
%ax.YAxis.FontSize = 4;
ylim([0,25])
yticks('manual')
yticks(0:1:25)
yticklabels(rightlabelstr')

% for thr = 1:size(U_s_log,1)
%     for k = 1:10:t_seconds
%         yline(thr, 'LineWidth', 0.75, 'Color','b' )
%         if U_s_log(thr,k) == 1
%             plot(k,thr,'r.','LineWidth',0.75)
%         end
%     end
% %plot(1:t_seconds, U_s_log,'LineWidth',0.50)
% end
title('Thruster Activity','FontSize',12)
xlabel('Time (seconds)','FontSize',14)

figure
plot(1:t_seconds,b2rquat_log,'LineWidth',1)
title('Body to Reference Quaternion vs. Time','FontSize',12)
xlabel('Time (seconds)','FontSize',14)
ylabel('Quaternion','FontSize',14)
legend('q0','q1','q2','q3','FontSize',12)

figure 
plot(1:t_seconds,deltarlog,'LineWidth',1)
title('Difference between Reference and Spacecraft Position','FontSize',12)
xlabel('Time (seconds)','FontSize',14)
ylabel('Distance (meters)','FontSize',14)
legend('x','y','z','FontSize',12)

figure 
plot(1:t_seconds,F_command_log,'LineWidth',1)
title('Commanded Translational Force','FontSize',12)
xlabel('Time (seconds)','FontSize',14)
ylabel('Force (newtons)','FontSize',14)
legend('x','y','z','FontSize',12)

figure
plot3(r_des_log(:,1),r_des_log(:,2),r_des_log(:,3),'og',r_dep_log(:,1),r_dep_log(:,2),r_dep_log(:,3),'--r','MarkerSize',0.5)
title('Deputy and Desired Deputy Trajectories')
legend('Desired', 'Deputy')

figure
plot(1:t_seconds,kep_log(:,1),'LineWidth',1)
title('Semimajor Axis vs. Time','FontSize',12)
xlabel('Time (seconds)','FontSize',14)
ylabel('Semimajor Axis (meters)','FontSize',14)

figure
plot(1:t_seconds,kep_log(:,2),'LineWidth',1)
title('Eccentricity vs. Time','FontSize',12)
xlabel('Time (seconds)','FontSize',14)
ylabel('Eccentricity','FontSize',14)

figure
plot(1:t_seconds,kep_log(:,3),'LineWidth',1)
title('Inclination vs. Time','FontSize',12)
xlabel('Time (seconds)','FontSize',14)
ylabel('Inclination (degrees)','FontSize',14)