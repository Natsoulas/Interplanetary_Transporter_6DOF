function res = attitude_dynamics_modelode5(t,z0,IC)
%attitude dynamics numerical integration MOI T_control
%state variable shall be Z = [quaternion,angularvel] dZ = [quaterniondot,angular accel]
%Initial state variable is Z0
% Inputs:
%   z0 (1x4 float): Initial conditions in order: quaternion (q0,q1,q2,q3),
%   angular velocity w (w1,w2,w3)
% Outputs:
%   t (1000x1 float): Array of times
%   res (1000x7 float): Array of state values corresponding to each time

res = ode5(@attitude_eom,linspace(0.001,t.span,t.steps),z0,IC);

    function dz = attitude_eom(t,z,IC)
        %z = [q0,q1,q2,q3,w1,w2,w3];
        %mrp = quat2mrp([z(1),z(2),z(3),z(4)]);
        q = [z(1);z(2);z(3);z(4)];
        q = q/norm(q);
        w = [z(5);z(6);z(7)];
        wdot = IC.MOI\(-[0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0]*IC.MOI*w + IC.T_control(:,round(1000*t)));
        wdot1 = wdot(1);
        wdot2 = wdot(2);
        wdot3 = wdot(3);
        qdot = 0.5*[0,-w(1),-w(2),-w(3);w(1),0,w(3),-w(2);w(2),-w(3),0,w(1);w(3),w(2),-w(1),0]*q;
        qdot0 = qdot(1);
        qdot1 = qdot(2);
        qdot2 = qdot(3);
        qdot3 = qdot(4);
        dz = [qdot0;qdot1;qdot2;qdot3;wdot1;wdot2;wdot3];
    end

end