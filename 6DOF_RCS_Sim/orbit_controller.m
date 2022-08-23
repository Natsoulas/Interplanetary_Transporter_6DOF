function [F_command,Deltar,Deltar_d,A,B] = orbit_controller(z,IC)
%z = [q0,q1,q2,q3,w1,w2,w3,rx,ry,rz,vx,vy,vz,rxd,ryd,rzd,vxd,vyd,vzd,rxc,ryc,rzc,vxc,vyc,vzc];
r_deputy = z(8:10)'- Earthrfunc(IC.t_orb,IC);
v_deputy = z(11:13)' - Earthvfunc(IC.t_orb,IC);
r_deputy_desired = z(14:16)'- Earthrfunc(IC.t_orb,IC);
v_deputy_desired = z(17:19)'- Earthvfunc(IC.t_orb,IC);
% r_chiefmag = norm(r_chief);
% r_depmag = norm(r_deputy);
% rdepdesiredmag = norm(r_deputy_desired);
f_kep = @(r,mu) -(mu/(norm(r)^3))*r ;
Deltar = r_deputy - r_deputy_desired;
Deltar_d = v_deputy - v_deputy_desired;
A = f_kep(r_deputy,IC.mu_auto);
B = f_kep(r_deputy_desired,IC.mu_auto);
F_com = -(A-B) - IC.K1*Deltar - IC.K2*Deltar_d;
F_command = F_com*IC.m_sc;
end