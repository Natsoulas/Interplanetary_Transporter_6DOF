%Control_Loop.m copy for ode45 implementation. Discontinous torque vector
%is splined into a continuous vector so that a varying timestep integrator
%(ode45) can quickly solve it
clc
close all
clear all
%Define IC's
%(junkbox)==============================================================================================
%%%Describe and propagate orbit of celestial bodies to be globally
IC.torcapx = 16000;
IC.torcapy = 16000;
IC.torcapz = 3500;
IC.fcap = 800;
%%%available
IC.m_sc = 90000;
IC.mu = 1.32712440042E20; %sun Gm in m^3/day^2 %gravitational parameter of sun
% [Afun, Bfun, Mfun, rfun, vfun, rfunv] = KeplerFun();
%  % Define Planet orbits with Keplerian elements
 a_e = (1.495468303960526E+08)*1000; %km 2 meters
 a_moon = (3.810068005390199E+05)*1000;%meters, [MAKE CENTER EARTH]
 a_mars = (2.279419507091813E+08)*1000;% km 2 meters
 e_e = 1.699803590732099E-02; % 
 e_moon = 5.875998505883045E-02;%average
 e_mars = 9.358077935297598E-02;
 I_e = 4.787963974336582E-03*pi/180; %rad
 I_moon = deg2rad(5.273587333467154E+00);%from ecliptic
 I_mars = deg2rad(1.846477690234122E+00); %degrees input
 w_e = 2.700238217175390E+02*pi/180;
 w_moon = deg2rad(1.909718909133165E+02);%
 w_mars = deg2rad(2.867965399000748E+02);%
 O_e = 1.917131048127610E+02*pi/180;
 O_moon = deg2rad(8.143902414246456E+01);%
 O_mars = deg2rad(4.943608223489434E+01);
 t_p_e = 2465791.336754327174; %JD5
 t_p_moon = 2465974.859928430058;%same situation as below
 t_p_mars = 2465934.922529480420;%in JD number figure out how to convert to seconds
 IC.mu_e = 398600.435436*(1000^3); %km^3/s^2 convert to m3/s2
 IC.mu_moon = 4902.800066*(1000^3);%km^3/s^2 converted to m3/s2
 IC.mu_mars = 42828.375214*(1000^3);%km^3/s^2 convert to m3/s

% %graph of orbits
 
% r_es = rfunv(a_e,e_e,I_e,w_e,O_e,...
%    invertKepler(Mfun(n_e,0:10000000:2*pi/n_e,t_p_e),e_e));
%  r_moon = rfunv(a_moon,e_moon,I_moon,w_moon,O_moon,...
%    invertKepler(Mfun(n_moon,0:10000:2*pi/n_moon,t_p_moon),e_moon));
%  r_mars = rfunv(a_mars,e_mars,I_mars,w_mars,O_mars,...
%    invertKepler(Mfun(n_mars,0:100000000:2*pi/n_mars+100,t_p_mars),e_mars));

%  %plot results
%  figure(1)
%  clf
%  plot3(r_es(1,:),r_es(2,:),r_es(3,:),r_moon(1,:),r_moon(2,:),r_moon(3,:),...
%    r_mars(1,:),r_mars(2,:),r_mars(3,:), 'LineWidth',2)
%     axis equal
%%
%Main engines burn data: Choose body frame force and burn time/duration
IC.t_burnstart = 0;
IC.t_burnend = 60*30;
IC.burnforce = [0;0;0];
%%
 IC.r_soi_e = a_e*((IC.mu_e)/IC.mu)^(2/5);
 IC.r_soi_moon = a_moon*((IC.mu_moon)/IC.mu)^(2/5);
 IC.r_soi_mars = a_mars*((IC.mu_mars)/IC.mu)^(2/5);
%%%Describe spacecraft geometry/inertia/thrusters
IC.MOI = diag([1.34E5,7.6E4,2.82E4]); % moment of inertia matrix for the spacecraft
%IC.thrusterdata = load('thruster_plac_data.mat');
IC.thrusterdata = load('sc_thrusterdatagood.mat');
IC.thrusterdata = IC.thrusterdata.MTAS_thrusterdata;
IC.R = 4.0;
IC.maxthrusterson = 5;
%IC.thrusterdata = IC.thrusterdata.Thruster_PlacementData;
%IC.thrusterdata = [zeros(48,1),IC.thrusterdata(:,2:end)];
%IC.thrusterdata(:,8) = 25;
%IC.thrusterdata = str2double(IC.thrusterdata);
IC.CG = [0;0;5.5]; %center of gravity/center of mass of spacecraft
%%%Rotational dynamics parameters
IC.CNTRL_MODE = 'slew'; %'translate'

% IC.N_C_B = rotx(45); % DCM for body to inertial: identity matrix to say they are the same initially
% IC.N_C_R = rotx(45); % DCM for ref to inertial: 90 degree rotation matrix about the z-axis

%IC.w_bn0 = transpose([0 0 0]);
%IC.qb0_N = transpose([1 0 0 0]);
%%%reference------------------------------------------
%ref.w_rn_R = [0;0;0]; %reference angular elocity of r/n in R frame
%ref.wdotrn_R = [0;0;0]; %time derivative of element above
%%%----------------------------------------------------
%delta---updates throughout loop
%IC.delta_w_bn0 = IC.w_bn0 - inv(IC.N_C_B)*IC.N_C_R*ref.w_rn_R;
%external torques (gravity, drag, etc.)
IC.L_external = 0; %external torques in EOM (0 for now)
% mass flow for a thruster (assuming all are sized the same)
IC.mdot = 200/(295*9.8); %kg/s
%gains for controllers and pwpf modulator
IC.K1 = 15/2500000;%0.0000032*eye(3); %was 3,2
IC.K2 = 22/2500000;%0.0000044*eye(3);%was 4.5,3.5
IC.K = 5000; %check boulder slides %was 10000
IC.P_matrix = IC.K*400*[1,0,0;0,4,0;0,0,4]; %3rd element of diagonal needs to be much larger for com z offset
IC.C = 1; %pwpf command signal (set at 1 to accept burn durations from simplex)
IC.K_p = 3; % proportional tuning gain for pwpf
IC.K_m = 2.5; % tuning gain for pwpf
IC.T_m = 0.28; %tuning gain for pwpf
IC.U_on = 0.9;%0.9; %tuning parameter for pwpf (schmitt trigger)
IC.U_off = IC.U_on/7; %tuning parameter for pwpf (schmitt trigger)
%tolerance for Schmitt Trigger and time structure for numint
IC.tol = 0.004; %tolerance
t.span = 3600;
t.steps = 3600;
t.microsteps = 1000;
t.microspan = 1;
%%% Initial state
mu_init = IC.mu_e; %meters and seconds please (initial grav param for scenario) %central body: earth
%Ra = 399000*1000; %meters
JD_init = juliandate(datetime('01-Jul-2039'));
%Rp = 6778*1000;
% a_init = 42164000; %meters
%  e_init = 0;
%  i_init = 0;%degrees
%  o_init = 0;%
%  O_init = 0;%
%  nu_init =40;% degrees
%  truelon = 20;
%  [r_chief,v_chief] = keplerian2ijk(a_init,e_init,i_init,O_init,o_init,nu_init,'truelon',truelon); %enter orbital elemetns for situation's desired orbit
%  r_chief = rotx(-23.4)*r_chief + Earthr(JD_init);
%  v_chief = rotx(-23.4)*v_chief + Earthv(JD_init)/86400;
 r_chief = [-6809.0233*1000;-18634.711*1000;-1650.2238*1000];
v_chief = [5.13918442*1000;3.40251526*1000;0.30131465*1000];
[a_i,e_i,i_i,O_i,o_i,nu_i] = ijk2keplerian(r_chief,v_chief);
[r_peri,v_peri] = keplerian2ijk(a_i,e_i,i_i,O_i,o_i,0); %inertial position and velocity at perigee

IC.r_peri = r_peri + Earthrfunc(JD_init,IC);
IC.v_peri = v_peri + Earthvfunc(JD_init,IC);

nu_current = 330;

[r_chiefy,v_chiefy] = keplerian2ijk(a_i,e_i,i_i,O_i,o_i,nu_current); 



r_chief = r_chiefy + Earthrfunc(JD_init,IC);
v_chief = v_chiefy + Earthvfunc(JD_init,IC);

% r_rot1 = vrrotvec(r_chief/norm(r_chief),IC.r_peri/norm(IC.r_peri));
% r_c2p = vrrotvec2mat(r_rot1);
% r_rot2 = vrrotvec(v_chief/norm(v_chief),IC.v_peri/norm(IC.v_peri));
% r_c2p2 = vrrotvec2mat(r_rot2);
% rotsum = r_c2p2*r_c2p;
IC.t_orb = JD_init;
radial = r_chief - Earthrfunc(IC.t_orb,IC);
rhat = radial/norm(radial);
   
velly = v_chief-Earthvfunc(IC.t_orb,IC);
vhat = (velly/norm(velly));
h = cross(radial,velly);
hhat = h/norm(h);

thehat = -cross(vhat,hhat);

IC.N_C_B = [thehat,hhat,vhat];

%disp(IC.N_C_B)
IC.N_C_R = [thehat,hhat,vhat];
CHEK = IC.N_C_R'*IC.N_C_B;
CHEK2 = dcm2quat(CHEK);
IC.qr0_B = dcm2quat(IC.N_C_B);
r_des = r_chief; 
v_des = v_chief;
r_sc = r_chief+[0;0;0]; 
v_sc = v_chief;%sqrt(mu_init*((2/norm(r_sc)) - 1/a_init))*(v_des/norm(v_des)); %vis-viva equation
nud = norm(h)/(norm(radial)^2);
worb0 = ((1+e_i*cosd(nu_i))/(1 + e_i^2 + 2*e_i*cosd(nu_i)))*nud*hhat;
w_0 =inv(IC.N_C_B)*worb0;
%z = [q0,q1,q2,q3,w1,w2,w3,rx,ry,rz,vx,vy,vz,rxd,ryd,rzd,vxd,vyd,vzd,nu,nudot];
z0 = [IC.qr0_B(1), IC.qr0_B(2),IC.qr0_B(3),IC.qr0_B(4),w_0(1),w_0(2),w_0(3),r_sc(1),r_sc(2),r_sc(3),v_sc(1),v_sc(2),v_sc(3),r_des(1),r_des(2),r_des(3),v_des(1),v_des(2),v_des(3)];
z = z0; %intializes z as z0

IC.mu_auto = IC.mu_e;
%log z (state)
z_log = [];
r_dep_log = [];
r_des_log = [];
kep_log = [];
%%%reformatting thruster data (ignore pls)
thrusterin4simplex = [];
unsplinedT_log = [];
unsplinedF_log = [];
propellant_used = 0; %kilograms
for k = 1:1:size(IC.thrusterdata,1)
        angle = IC.thrusterdata(k,3);
        [X,Y] = azim2cartbody(angle,IC.CG,IC.R);
        h = IC.thrusterdata(k,4);
        ZwrtCOM = h - IC.CG(3);
        Z = ZwrtCOM;
        U = IC.thrusterdata(k,5);
        V = IC.thrusterdata(k,6);
        W = IC.thrusterdata(k,7);
        force = IC.thrusterdata(k,8);
        thrusterin4simplex = [thrusterin4simplex; k,X,Y,Z,U,V,W,force];
end
%ref.Force = [0;0;0]; %no translational force for now please
%======================================================================================================
%establish control loop
t_seconds = 0;
B2Rquat = [3; 2; 1; 4];
Torque_res_log = [];
Force_res_log = [];
controllerout_log = [];
U_s_log = [];
sollog = [];
b2rquat_log = [];
F_command_log =[];
deltarlog = [];
deltar = [0;0;0];
deltar_d = [0;0;70];
A_log = [];
B_log = [];
%norm(B2Rquat(2:4)) > 0.04 && 
IC.t_seconds = 0;
while round(nu_current) ~=  30 %&& t_seconds < 10000%norm(deltar) > 5 && norm(deltar_d) > 3
    IC.t_orb = t_seconds/86400 + JD_init;
    if mod(t_seconds,100) == 0
        disp('Seconds Passed: ')
        disp(t_seconds)
    end
    %starts with parameters fed into attitude controller
    %z = z0 for sake of input
    %[L_command,sig,B2Rquat] = attitude_controller(z,IC,ref);
    [L_command,sig,B2Rquat] = attitude_controller_rntref(z,IC);
    L_command = sign(L_command).*[min(abs(L_command(1)),IC.torcapx),min(abs(L_command(2)),IC.torcapy),min(abs(L_command(3)),IC.torcapz)]';
    %[r_Lcom,c_Lcom] = size(L_command);
    %L_command = zeros(r_Lcom,c_Lcom);
    [F_command,deltar,deltar_d,A_a,B_b] = orbit_controller(z,IC);
    capx = min(F_command(1),IC.fcap);
    F_command = sign(F_command).*[min(abs(F_command(1)),IC.fcap),min(abs(F_command(2)),IC.fcap),min(abs(F_command(3)),IC.fcap)]';
    A_log = [A_log, A_a];
    B_log = [B_log, B_b];
    if mod(t_seconds,100) == 0
        disp('F_command: ')
        disp(F_command')
        disp('L_command: ')
        disp(L_command')
    end
    F_command_log = [F_command_log,F_command];
    deltarlog = [deltarlog,deltar];
    b2rquat_log = [b2rquat_log; B2Rquat];
    %disp(L_command)
    if mod(t_seconds,100) == 0
        disp('B2R Quat: ')
        disp(B2Rquat)
    end
    controllerout_log = [controllerout_log,L_command];
    % simplex
    B_C_N = quat2dcm(z(1:4));
    solution = simplexrunner(thrusterin4simplex,L_command,B_C_N'*F_command,IC.CG,IC.CNTRL_MODE);
    sollog = [sollog, solution(:,1)];
    % PWPF
    U_s = zeros(size(solution,1),1000);
    for yi = 1:1:size(solution,1)
    t_s = solution(yi,1);
    t_s = round(t_s*1000)/1000;%rounds seconds to nearest millisecond
    if t_s > 1.00
        t_s = 1.00;
    end
    microsteps = round(t_s/0.001); %tells number of milliseconds which makes the stepsize 1 millisecond
    if t_s > 0
        [u,DC,f_o] = PWPF_Run(IC.C,IC.K_p,IC.K_m,IC.T_m,IC.U_on,IC.U_off,t_s,microsteps,IC.tol);
        for ji = 1:size(u,2)
            U_s(yi,ji) = u(ji);
        end
    end
%     disp('DC')
%     disp(DC)
%     disp('f_o')
%     disp(f_o)
    end
    U_s_log = [U_s_log, U_s(:,1:end)];
    % PWPF output is then processed into torques given thruster data
    %%bite of foreign code
        %compute total torque vectors for all unique sets of n-thrusters
    dimlong = size(IC.thrusterdata,1);
    thrusterno = 1:dimlong;
    Torque_ms = zeros(3,1000);
    for millisec = 1:1000
        %loop through which millisecond you are on.
        Tmilli = [0;0;0];
        Fmilli = [0;0;0];
        for k = 1:dimlong
            angle = IC.thrusterdata(k,3);
            [X,Y] = azim2cartbody(angle,IC.CG,IC.R);
            h = IC.thrusterdata(k,4);
            ZwrtCOM = h - IC.CG(3);
            Z = ZwrtCOM;
            r = [X;Y;Z];
            U = IC.thrusterdata(k,5);
            V = IC.thrusterdata(k,6);
            W = IC.thrusterdata(k,7);
            force = IC.thrusterdata(k,8);
            F = force*[U;V;W];
            if U_s(k,millisec) == 1
                Tconcat = cross(r,F);
                Fconcat = F;
                propellant_used = propellant_used + 0.001*IC.mdot;
            else
                Tconcat = [0;0;0];
                Fconcat = [0;0;0];
            end
            %sums all the produced torques for the particular millisecond
            %but only if the partifular thruster is on.
            Tmilli = Tmilli + Tconcat;
            Fmilli = Fmilli + Fconcat;
        end
        Torque_ms(:,millisec) = Tmilli;
        Force_ms(:,millisec) = Fmilli;
    end
    %%
    %spline force array to be continuous
    unsplinedF_log = [unsplinedF_log,Force_ms(:,1:10:end)];
    x_spl = 0.001:0.001:1;
    y_spl_1 = Force_ms(1,:);
    xx_spl = 0.001:0.001:1;
    yy_spl_1 = spline(x_spl,y_spl_1,xx_spl);
    y_spl_2 = Force_ms(2,:);
    yy_spl_2 = spline(x_spl,y_spl_2,xx_spl);
    y_spl_3 = Force_ms(3,:);
    yy_spl_3 = spline(x_spl,y_spl_3,xx_spl);
    Force_ms = [yy_spl_1;yy_spl_2;yy_spl_3];
    Force_ms_log = Force_ms(:,1:100:end);
    %spline torque array to be continuous
    unsplinedT_log = [unsplinedT_log,Torque_ms(:,1:10:end)];
    x_spl = 0.001:0.001:1;
    y_spl_1 = Torque_ms(1,:);
    xx_spl = 0.001:0.0001:1;
    yy_spl_1 = spline(x_spl,y_spl_1,xx_spl);
    y_spl_2 = Torque_ms(2,:);
    yy_spl_2 = spline(x_spl,y_spl_2,xx_spl);
    y_spl_3 = Torque_ms(3,:);
    yy_spl_3 = spline(x_spl,y_spl_3,xx_spl);
    Torque_ms = [yy_spl_1;yy_spl_2;yy_spl_3];
    Torque_ms_log = Torque_ms(:,1:100:end);
    % torques summed and fed into ode5
    IC.T_control = Torque_ms;
    IC.F_control = Force_ms; %the edit
    tspan = 1;
    tsteps = 1000;
    [t,res] = orbit_attitude_numint_VEL_FIRE(tspan,tsteps,z,IC);
    z_update_2controller = res(end,:);
    %disp(z_update_2controller(1,5:7))
    % ode5 sends state back through loop for controller to deal with
        %updates state variables in state row vector z
    z_log = [z_log; res(end,:)];
    
    z = z_update_2controller;
    r_dep_log = [r_dep_log;z_update_2controller(1,8:10)-Earthrfunc(IC.t_orb,IC)'];
    r_des_log = [r_des_log;z_update_2controller(1,14:16)-Earthrfunc(IC.t_orb,IC)'];
    %t.span = 3600;
    %t.steps = 3600;
    Torque_res_log = [Torque_res_log, Torque_ms_log];
    Force_res_log = [Force_res_log, Force_ms_log];
   
    r_current = z(8:10)'-Earthrfunc(IC.t_orb,IC);
    v_current = z(11:13)' - Earthvfunc(IC.t_orb,IC);
    orbelarray = coord_set_convert('cart', [r_current',v_current'], 0.0001, IC.mu_e);
    nu_current = rad2deg(orbelarray(6));
    if mod(t_seconds,100) == 0
        disp('True Anomaly: ')
        disp(rad2deg(orbelarray(6)))
    end

    orbelarray(3:6) = rad2deg(orbelarray(3:6));
    kep_log = [kep_log; orbelarray];

     t_seconds = t_seconds + 1;
    IC.t_seconds = t_seconds;
end