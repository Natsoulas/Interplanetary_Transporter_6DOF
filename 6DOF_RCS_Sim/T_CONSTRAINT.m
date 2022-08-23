function [F,G] = T_CONSTRAINT(X,f_args)
% This is a user defined function that spits out the current solution cost
% (F) and constraint vector (G) that are used in the main simps function to
% determine if a solution sucks


% f_args = {A,b,ISP,G_mat,BND_PEN,PROP_PEN,TNUM_PEN,TNUM_LIM,X_LB,X_UB,say_stuff,ALT_PEN,alt_A,alt_b,G_mat_alt};
say_stuff = f_args{11};


X = X';
A = f_args{1};
b = f_args{2};
A_alt = f_args{13};
b_alt = f_args{14};
ISP = f_args{3};
W_G = f_args{4};
W_G_alt = f_args{15};

F = (X')*ISP; % NOTE: This cost is augmented with more crap in the main 'simps'
%               function (starting on line 302 in simps). This is the
%               baseline cost of propellant usage
G = W_G*(A*(X) - b); % This is the cost of not achieving the desired force/torque
G_alt = W_G_alt*(A_alt*(X)-b_alt);
%                      scaled by W_G (which is G_mat in the main script)

% This will output a bunch of information to the command window each time
% this function is called, which is many times, so I recommend leaving
% 'say_stuff' to false to avoid that
if say_stuff

    % The code in here effectively replicates what's going on in the
    % augmented cost portion of 'simps' (starting line 302 in simps)
    BND_PEN = f_args{5};
    PROP_PEN = f_args{6};
    TNUM_PEN = f_args{7};
    TNUM_LIM = f_args{8};
    %ALT_PEN = f_args{12};
    X_LB = f_args{9};
    X_UB = f_args{10};

    LBND_CONSTR = X_LB - X;
    UBND_CONSTR = X - X_UB;

    COUNT_PEN = 1;
    if numel(X(X > 0.00001)) > TNUM_LIM %changed x to X thurs 7/28
        COUNT_PEN = TNUM_PEN*(numel(X(X > 0.00001))-TNUM_LIM);
    end
    % Should heavily penalize cost if using x outside provided bounds or
    % violating equality constraints
    F1 = PROP_PEN*F; % Prop usage penalty
    %F2 = 1e12*sum(LBND_CONSTR(LBND_CONSTR > -1e-10)); % Lower bound violation penalty
    F2 = 1e8*sum(abs(LBND_CONSTR(LBND_CONSTR > -1e-10)));
    F3 = BND_PEN*sum(UBND_CONSTR(UBND_CONSTR > 0)); % Upper bound penalty
    F4 = sum(abs(G)); % Target force and torque penalty
    F5 = COUNT_PEN;%*F; % Thruster count penalty
    F6 = sum(abs(G_alt)); %target of alternate of mode penalty (if in slew, penalizes for missing trans force) etc.
    %F_tot = F1 + F2 + F3 + F4 + F5 + F6; % Target force and torque penalty

    fprintf('Prop Cost: %1.5e\n',F1)
    fprintf('LB Cost  : %1.5e\n',F2)
    fprintf('UB Cost  : %1.5e\n',F3)
    fprintf('TRGT Cost: %1.5e\n',F4)
    fprintf('NUM Cost : %1.5e\n',F5)
    fprintf('ALT Cost: %1.5e\n',F6)
    % NOTE: If you want to add or reduce a cost term, then you have to
    % change it here AND over in simps
end

G = G';
end