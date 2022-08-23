function [L_command,sig,quaty] = attitude_controller(z,IC)
%takes pwpf modulator gains as input as well as state z for each timestep
%to calculate the control output torque of L_command so that it can be
%propagated through the control loop.
%keep in mind that "w" is angular velocity in the B (body) frame.
%keep in mind that z = [q0 q1 q2 q3 w1 w2 w3];
quat = z(1:4);
IC.N_C_B = quat2dcm(quat)'; %I flipped it cuz quat2dcm outputs B_C_N %put back inv if error

r_des =z(14:16)';
radial = r_des - Earthrfunc(IC.t_orb,IC);
rhat = radial/norm(radial);
   
velly = z(17:19)'-Earthvfunc(IC.t_orb,IC);
vhat = (velly/norm(velly));
h = cross(radial,velly);
hhat = h/norm(h);

thehat = cross(vhat,hhat);

obelvec = coord_set_convert('cart', [radial',velly'], 0.0001, IC.mu_e);

e_c = obelvec(2);
nu_c = rad2deg(obelvec(6));
        nud = norm(h)/(norm(radial)^2);
        nudd = -2*(norm(velly)/norm(radial))*nud;
        orbw = ((1+e_c*cosd(nu_c))/(1 + e_c^2 + 2*e_c*cosd(nu_c)))*nud*hhat;
        orbwd = (((1+e_c*cosd(nu_c))/(1 + e_c^2 + 2*e_c*cosd(nu_c)))*nudd - ((e_c*(1-e_c^2)*sind(nu_c))/((1 + e_c^2 + 2*e_c*cosd(nu_c))^2)*nud*nud))*hhat;

IC.N_C_R = [thehat,vhat,hhat]';
w = transpose(z(5:7));
wcross = [0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0];
quaty = dcm2quat(inv(IC.N_C_R)*IC.N_C_B);
sig = quat2mrp(transpose(quaty)); %symbol for MRP
%



IC.w_bn0 = w;
%get rid of IC on N_C_B
w_rn_R = orbw;
w_ref_dot = orbwd;
delta_w= IC.w_bn0 - inv(IC.N_C_B)*IC.N_C_R*w_rn_R; %difference between ref and currrent angular velocities - expression from cu boulder slides
L_command = -IC.K*sig - IC.P_matrix*delta_w + IC.MOI*(inv(IC.N_C_B)*IC.N_C_R*w_ref_dot - wcross*inv(IC.N_C_B)*IC.N_C_R*w_rn_R) + wcross*IC.MOI*w - IC.L_external; %expression from cu boulder slides
end