function [t,res] = orbit_attitude_numint_VEL_FIRE(tspan,tsteps,z0,IC)
%attitude dynamics numerical integration
%state variable shall be Z = [quaternion,angularvel] dZ = [quaterniondot,angular accel]
%Initial state variable is Z0
% Inputs:
%   z0 (1x4 float): Initial conditions in order: quaternion (q0,q1,q2,q3),
%   angular velocity w (w1,w2,w3)
% Outputs:
%   t (1000x1 float): Array of times
%   res (1000x7 float): Array of state values corresponding to each time
%indices = 1:1:size(IC.T_control,2);

%anonymous functions for cowell's method numerical integration
% qf = @(r21,rj1) dot(r21,r21-2*rj1)/dot(rj1,rj1);
% fqf = @(q) q*(3+3*q+q^2)/(1 + (1+q)^(3/2));
%
% ode45 integration call
ode_options = odeset('RelTol',1e-10,'AbsTol',1e-11*ones(size(z0,1)));
[t,res] = ode45(@orbit_attitude_eom,linspace(0.001,tspan,tsteps),z0,ode_options);

    function dz = orbit_attitude_eom(t,z)

        %z = [q0,q1,q2,q3,w1,w2,w3,rx,ry,rz,vx,vy,vz,rxd,ryd,rzd,vxd,vyd,vzd,rxc,ryc,rzc,vxc,vyc,vzc];
        %rxd and rxc have d and c at the end meaning deputy and chief
        %respectively
        %mrp = quat2mrp([z(1),z(2),z(3),z(4)]);
        %orbit
        r21 = z(8:10); %spacecraft to sun position vector
        r21_d = z(14:16);
        r21mag = norm(r21);
        r31 = Earthrfunc(IC.t_orb+t/86400,IC); %earth to sun position vector
        r23 = r21-r31;
        r23_d = r21_d-r31;
        r23mag = norm(r23);
        r23dmag = norm(r23_d);

        r41 = Moonrfunc(IC.t_orb+t/86400,IC); %moon to sun position vector (vector additon of earth2sun and moon2eartn positons)
        r24 = r21-r41;
        r24mag = norm(r24);
        r51 = Marsrfunc(IC.t_orb+t/86400,IC); %mars to sun position vector
        r25 = r21-r51;
        r25mag = norm(r25);
%         qone = qf(r21, r31);
%         qone_d = qf(r21_d,r31);
%         qtwo = qf(r21, r41);
%         qthree = qf(r21, r51);
%         fq1 = fqf(qone);
%         fq2 = fqf(qtwo);
%         fq3 = fqf(qthree);
%         fq1_d = fqf(qone_d);
        Earthgrav = -(IC.mu_e)/r23mag^3*(r23);% + fq1*r31);
        Earthgrav_ad = -(IC.mu_e)/r23dmag^3*(r23_d);% + fq1_d*r31);
        Lunargrav = -(IC.mu_moon)/r24mag^3*(r24);% + fq2*r41);
        Martiangrav = -(IC.mu_mars)/r25mag^3*(r25);% + fq3*r51);
        %P_SRP = 0; %acceleration due to Solar Radiation Pressure
        %%
        %conditionals for toggles
        if r23mag <= IC.r_soi_e %|| t < 2.462847082058206e+06
            earth_toggle = 1;
        else
            earth_toggle = 0;
        end
        if r24mag <= IC.r_soi_moon
            moon_toggle = 1;
        else
            moon_toggle = 0;
        end
        if r25mag <= IC.r_soi_mars
            mars_toggle = 1;
        else
            mars_toggle = 0;
        end
        %%
        %z = [q0,q1,q2,q3,w1,w2,w3,rx,ry,rz,vx,vy,vz,rxd,ryd,rzd,vxd,vyd,vzd];
        vx = z(11);
        vy = z(12);
        vz = z(13);
        
%         vxc = z(23);
%         vyc = z(24);
%         vzc = z(25);
        
        vxd = z(17);
        vyd = z(18);
        vzd = z(19);
        Fcontrol = IC.F_control(:,round(t*size(IC.F_control,2)))/IC.m_sc;
        
        q = [z(1);z(2);z(3);z(4)];
        q = q/norm(q);
        B_C_N = quat2dcm(q');
        BURN = timedburn(IC.t_burnstart,IC.t_burnend,IC.burnforce,IC.t_seconds + t);
        a = -(IC.mu/r21mag^3)*r21 + earth_toggle*Earthgrav + moon_toggle*Lunargrav + mars_toggle*Martiangrav + Fcontrol + B_C_N'*BURN/IC.m_sc;
        ax = a(1);
        ay = a(2);
        az = a(3);
        
       
        %disp(a')

        %translational chief deputy
%         r_chief = z(20:22)';
        %r_dep = z(8:10)';% r_dep_d = z(14:16)';
        

        %ac = -(IC.mu/(norm(r_chief)^3))*r_chief;

        ad = -(IC.mu/(norm(r21_d)^3))*r21_d + earth_toggle*Earthgrav_ad + ([vxd;vyd;vzd]/norm([vxd;vyd;vzd]))*norm(BURN)/IC.m_sc;
        

  
%         axc = ac(1);
%         ayc = ac(2);
%         azc = ac(3);
%         if mod(IC.t_seconds,10)
%             disp(a-ad)
%         end
        axd = ad(1);
        ayd = ad(2);
        azd = ad(3);
        %attitude
%         q = [z(1);z(2);z(3);z(4)];
%         q = q/norm(q);
         w = [z(5);z(6);z(7)];
        
        % gravity gradient torque
%         if r23mag > r24mag && %change logic
%             % earth defines nadir direction
%             nadir = r23/r23mag;
%         elseif moon_toggle == 1
%             %moon defines nadir direction
%             nadir = r24/r24mag;
%         elseif mars_toggle == 1
%             %mars defines nadir direction
%             nadir = r25/r25mag;
%         end
        possibleradials = [r23,r24,r25];
        [minmag,I] = min([r23mag,r24mag,r25mag]);
        radial = possibleradials(:,I)/minmag;
        %%%Copy below for orbwxydsdfa
%         rpos = possibleradials(:,I);
%         rhat = radial/norm(radial);
%    
%         velly = [[vx;vy;vz]-Earthvfunc(IC.t_orb+t/86400,IC), [vx;vy;vz]-Moonvfunc(IC.t_orb+t/86400,IC),[vx;vy;vz]-Marsvfunc(IC.t_orb+t/86400,IC) ];
%         vhat = (velly(:,I)/norm(velly(:,I)));
%         h = cross(rpos,velly(:,I));
%         hhat = h/norm(h);
%         
        %tangential direction
        %tangent = [vx;vy;vz]/norm([vx;vy;vz]);
        %normal = cross(radial,tangent);
        %inertial frame is heliocentric eccliptic whereas RNT/RSW is an
        %intermediate and can be defined in inertial
        nadir_N = -radial;
        %B_C_N = quat2dcm(q');
        nadir_B = B_C_N*nadir_N;
       
%         Xangle = N_eul_B(3);
%         Yangle = N_eul_B(2);
        R = minmag;
        mews = [IC.mu_e,IC.mu_moon,IC.mu_mars];
        mu = mews(I);
        IC.mu_auto = mu;
%         [a_c,e_c,i_c,O_c,o_c,nu_c] = ijk2keplerian(rpos,velly(:,I));
%         nud = norm(h)/(norm(rpos)^2);
%         nudd = -2*(norm(velly(:,I))/norm(rpos))*nud;
%         orbw = ((1+e_c*cosd(nu_c))/(1 + e_c^2 + 2*e_c*cosd(nu_c)))*nud*hhat;
%         orbwx = orbw(1);
%         orbwy = orbw(2);
%         orbwz = orbw(3);
%         orbwd = (((1+e_c*cosd(nu_c))/(1 + e_c^2 + 2*e_c*cosd(nu_c)))*nudd - ((e_c*(1-e_c^2)*sind(nu_c))/((1 + e_c^2 + 2*e_c*cosd(nu_c))^2)*nud*nud))*hhat;
%         orbwxd = orbwd(1);
%         orbwyd = orbwd(2);
%         orbwzd = orbwd(3);

        Tgg = cross((3*mu/(R^3))*nadir_B,IC.MOI*nadir_B);
        if mod(IC.t_seconds,100) == 0 && t == 0.001
            disp('Gravity gradient torque: ')
            disp(Tgg')
        end

        %turn off don't foergetr
        Tgg = [0;0;0];

        wdot = IC.MOI\(-[0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0]*IC.MOI*w + IC.T_control(:,round(t*size(IC.T_control,2)))+ Tgg);
        wdot1 = wdot(1);
        wdot2 = wdot(2);
        wdot3 = wdot(3);
        qdot = 0.5*[0,-w(1),-w(2),-w(3);w(1),0,w(3),-w(2);w(2),-w(3),0,w(1);w(3),w(2),-w(1),0]*q;
        qdot0 = qdot(1);
        qdot1 = qdot(2);
        qdot2 = qdot(3);
        qdot3 = qdot(4);
        
        
        
        dz = [qdot0;qdot1;qdot2;qdot3;wdot1;wdot2;wdot3;vx;vy;vz;ax;ay;az;vxd;vyd;vzd;axd;ayd;azd];
    end

end