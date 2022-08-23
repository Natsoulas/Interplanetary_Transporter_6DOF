% Simplex Method FUNCTION SCRIPT for modularity
function solution = simplexrunner(t_table,L_ref,F_ref,rCG_B,CNTRL_MODE)

if norm(L_ref) == 0 && norm(F_ref) == 0
    solution = zeros(24,1);
    return
end
curr_dir = pwd;

% Load in sample thruster table. Feel free to generate your own using the
% placement tool
%t_table = csvread([curr_dir,'/T_table.csv'],1,0);
N = size(t_table,1);

%==========================================================================
% User inputs


%L_ref = [1;1;0]; % Desired torque components in body frame
%F_ref = [0;0;0]; % Desired translational force components in body frame

%rCG_B = [0.1;...
%    0.25;...
%    0.00]; % Spacecraft CG referenced from bottom of spacecraft (so still
%            body frame, but origin is at bottom of spacecraft)

%==========================================================================

F = t_table(:,8);
rvec_B = t_table(:,2:4)'; % Thruster location in body frame
fvec_B = t_table(:,5:7)'; % Thruster force direction in body frame
Fvec_B = (t_table(:,5:7).*F)'; % Thruster force vector in body frame
Lvec_B = cross((rvec_B - rCG_B),Fvec_B); % thruster torque in body frame


%==========================================================================
% This portion looks at each thruster to determine which ones contribute to
% obtaining the desired torque and translational forces defined by L_ref
% and F_ref. The method considers the two goals separately, then uses a
% basic voting scheme to award "points" to thrusters that do one or both


% Projection of individual thruster force/torque in body frame directions
d_arr = [dot(cross((rvec_B - rCG_B),fvec_B),repmat([1;0;0],1,N));...
         dot(cross((rvec_B - rCG_B),fvec_B),repmat([0;1;0],1,N));...
         dot(cross((rvec_B - rCG_B),fvec_B),repmat([0;0;1],1,N))];
g_arr = [dot(fvec_B,repmat([1;0;0],1,N));...
         dot(fvec_B,repmat([0;1;0],1,N));...
         dot(fvec_B,repmat([0;0;1],1,N))];
Dx = d_arr(1,:); % Portion of thruster torque in x-direction
Dy = d_arr(2,:); % Portion of thruster torque in y-direction
Dz = d_arr(3,:); % Portion of thruster torque in x-direction
Gx = g_arr(1,:); % Portion of thruster force direction in x-direction
Gy = g_arr(2,:); % Portion of thruster force direction in y-direction
Gz = g_arr(3,:); % Portion of thruster force direction in z-direction

if strcmpi(CNTRL_MODE,'slew')
Fx = (Dx')*(Dx*(Dx'))^(-1)*dot(L_ref,[1;0;0]);
Fy = (Dy')*(Dy*(Dy'))^(-1)*dot(L_ref,[0;1;0]);
Fz = (Dz')*(Dz*(Dz'))^(-1)*dot(L_ref,[0;0;1]);

FFx = -ones(N,1);
FFy = FFx; FFz = FFx;

alt_A = Fvec_B;
alt_b = F_ref;

A = Lvec_B; % thruster forces and torques
b = L_ref; % desired force and torque
elseif strcmpi(CNTRL_MODE,'translate')
% Fx = (Dx')*(Dx*(Dx'))^(-1)*dot(L_ref,[1;0;0]);
% Fy = (Dy')*(Dy*(Dy'))^(-1)*dot(L_ref,[0;1;0]);
% Fz = (Dz')*(Dz*(Dz'))^(-1)*dot(L_ref,[0;0;1]);
Fx = -ones(N,1); Fy = Fx; Fz = Fx;
% This tells you what the force thruster N would need to have to provide
% the force necessary to get the x,y,z component of the desired force. So
% negative values imply the thruster should not be used (or should be
% flipped to point the other way...)
FFx = (Gx')*(Gx*(Gx'))^(-1)*dot(F_ref,[1;0;0]);
FFy = (Gy')*(Gy*(Gy'))^(-1)*dot(F_ref,[0;1;0]);
FFz = (Gz')*(Gz*(Gz'))^(-1)*dot(F_ref,[0;0;1]);

alt_A = Lvec_B;
alt_b = L_ref;

A = Fvec_B; % thruster forces and torques
b = F_ref; % desired force and torque
end



% These logical arrays indicate which thrusters provide a positive
% contribution to achieving the desired torque and force. A 1, or true,
% means the thruster helps by some extent
IIx = (Fx > 1); % x torque component
IIy = (Fy > 1); % y torque component
IIz = (Fz > 1); % z torque component
IIfx = (FFx > 0); % x force component
IIfy = (FFy > 0); % y force component
IIfz = (FFz > 0); % z force component

%==========================================================================



%%
TRGT_PEN = 1e6*ones(length(b),1); % Penalty for missing commanded force/torque
BND_PEN = 1e0; % Penalty for using more than upper bound for thruster
PROP_PEN = 1e0; % Penalty for prop usage
ALT_PEN = 1e3*zeros(length(alt_b),1);
TNUM_PEN = 1e5; % Penalty for using more than k-thrusters %niko edit - was 1e0 now is 1e1
TNUM_LIM = 8;%sum(F_ref > 0) + sum(L_ref > 0); % Maximum number of thrusters that 
%%                                             can be used before cost is
%                                             penalized
G_mat_alt = diag(ALT_PEN);
G_mat = diag(TRGT_PEN); % Penalty for missing desired torque/force applied evenly
%                         for all terms, so G_mat*(Ax-b)

X_LB = 0*ones(N,1); % Lower bound on thruster force. Strictly enforced
X_UB = ones(N,1); % Upper bound on thruster force. BND_PEN enforces this
XI = 1:N; % Acts as a switch. Each column with a '1' means the optimizer is 
%           free to select the optimal value to do the thing. If it is set
%           to '0', then that thruster will maintain its initial guess
%           value
ISP = ones(N,1); % This is the thruster mass flow rate. It is used in the
%                  cost function to determine a measure of prop usage

X_init = (X_UB - X_LB)*.05; % This is the initial guess. I assume all thrusters
% multiplier was .05 now is .01                            are initially on, but at a very low force
%                             value such that the optimizer is more likely
%                             to turn off ones that don't do much to help
%                             the cause while increasing those that do

say_stuff = false; % This is the 'verbose' flag that will output lots of info
%                    on the cost function terms and constraints

% This is an optional argument passed to simps. It provides simps and the
% constraint function with all of the parameters it needs to do the thing
f_args = {A,b,ISP,G_mat,BND_PEN,PROP_PEN,TNUM_PEN,TNUM_LIM,X_LB,X_UB,say_stuff,ALT_PEN,alt_A,alt_b,G_mat_alt};

% Here is a basic voting scheme for the thrusters. A thruster is awarded a
% point for each component of desired force/torque it helps to achieve
% based on the II arrays computed earlier.
% NOTE: A better way to do this might be to weigh it based on its specific
% contribution instead of a binary on/off. This info is held within the
% Fx/FFx, and so on, arrays
T_vote = sum([IIx,IIfx,IIy,IIfy,IIz,IIfz],2);
%T_vote = ones(24,1);
if sum(T_vote) == 0
    solution = zeros(24,1);
    return
end

XI(T_vote == 0) = []; % Don't use thrusters that don't help
ISP(T_vote > 1) = (1./(T_vote(T_vote > 1))).^2; % Reduce cost of using thrusters
%                                                 that help the most. NOTE:
%                                                 This could use some
%                                                 additional scaling so the
%                                                 optimizer actually cares.
%                                                 Maybe using the info in
%                                                 Fx/FFx and so on will
%                                                 help
X_init(T_vote == 0) = 0; % Set initial values of thrusters that don't help 
%                          to zero

options = zeros(1,14);
options(2) = 1e-4; % Tolerance
options(3) = 1e-4; % Tolerance
options(14) = 100*N; % Number of iterations
[x_opt,options] = simps(@T_CONSTRAINT,X_init,XI,options,X_LB,X_UB,f_args);
% [~,x_opt] = nma_simplex(A,b,ISP,false);

% After the run, turn 'verbose' on to take a close look at how well the
% constraints were satisfied as well as the solution vector and stuuuuff
say_stuff = false;
f_args = {A,b,ISP,G_mat,BND_PEN,PROP_PEN,TNUM_PEN,TNUM_LIM,X_LB,X_UB,say_stuff,ALT_PEN,alt_A,alt_b,G_mat_alt};

Fvec_LP = (t_table(:,5:7).*x_opt)'; % Compute what the force vec would be if the 
%                                     solution vec was used
Lvec_LP = cross((rvec_B - rCG_B),Fvec_LP); % Compute what the total torque vec
%                                            would be based on the solution
%                                            vec

solution = [x_opt,T_vote,ISP]; % [Solution vector, Thruster votes, Thruster prop cost]
%[fmin_test,~] = T_CONSTRAINT(x_opt',f_args); % Cost
A*(x_opt) - b; % Constraints
% disp('LOOK AT THIS THING')
% disp(transpose(A*(x_opt) - b))
% disp(';;;;;;;;;;;;;;;')
Fsum_LP = sum(Fvec_LP,2); % Force output (compare to desired F_ref)
Lsum_LP = sum(Lvec_LP,2);% Torque output (compare to desired L_ref)

% F(x_opt ~= 0) = x_opt(x_opt ~= 0);
% Fvec_B = (t_table(:,5:7).*F)'; % Thruster force vector in body frame
% Lvec_B = cross((rvec_B - rCG_B),Fvec_B); % thruster torque in body frame
% A = [Fvec_B;Lvec_B];
%
% [~,x_opt_test] = nma_simplex(A,b,f,false);
% x_opt_test
end
