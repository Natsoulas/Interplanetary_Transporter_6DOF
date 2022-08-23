function coord_set_out = coord_set_convert(coord_type, coord_set, dt, mu_val)
% Converts between Cartesian coordinates and Classical Orbital Elements for
% elliptical and hyperbolic orbits
% INPUTS:
%   - coord_type: String argument that specifies the type of coordinate set
%                 being supplied. Can be either 'cart' or 'orb el'
%   - coord_set : Set of elements to be converted to the other type. If
%                 orbital elements are supplied, the expected order is
%                 [SMA,ECC,INC,RAAN,ARGP,M0]. The angle values are expected
%                 to be in radians. The output state vectors will be in the
%                 units of the gravitational parameter provided.
%
%                 If a cartesian state vector provided, the corresponding
%                 orbital elements will be output

if strcmpi(coord_type,'cart') % Convert from cartesian to orb el
    
    r_vec = coord_set(1:3)';
    v_vec = coord_set(4:6)';
    h_vec = cross(r_vec,v_vec);
    
    r_mag = norm(r_vec);
    v_mag = norm(v_vec);
    h_mag = norm(h_vec);
    
    r_hat = r_vec/r_mag;
    e_vec = cross(v_vec,h_vec)/mu_val - r_hat;
    ECC = norm(e_vec);
    
    % Unit vectors in direction of parifocal coordinate frame axes
    h_hat = h_vec/h_mag;
    e_hat = e_vec/ECC;
    p_hat = cross(h_hat,e_hat);
    
    DCM_PN = [e_hat,p_hat,h_hat]'; % 3-1-3 Sequence Rotation Matrix
    
    % Use DCM to extract Euler angles (orbital elements in this case)
    % Limit angle ranges to [0 2pi]
    RAAN = wrapTo2Pi(atan2(DCM_PN(3,1),-DCM_PN(3,2)));
    INC = acos(DCM_PN(3,3));
    ARGP = wrapTo2Pi(atan2(DCM_PN(1,3),DCM_PN(2,3)));
    TA = atan2(dot(cross(e_hat,r_hat),h_hat),dot(e_hat,r_hat));
    SMA = mu_val*r_mag/(2*mu_val - v_mag^2*r_mag);
    
    E_tot = v_mag^2/2 - mu_val/r_mag; % Orbit type check
    if E_tot < 0 % Elliptical Case
        E = 2*atan(sqrt((1 - ECC)/(1 + ECC))*tan(TA/2));
        M = E - ECC*sin(E);
        M0 = M - sqrt(mu_val/(SMA^3))*dt;
        coord_set_out = [SMA,ECC,INC,RAAN,ARGP,wrapTo2Pi(M0)];
    elseif E_tot > 0 % Hyperbolic case
        H = 2*atanh(sqrt((ECC - 1)/(ECC + 1))*tan(TA/2));
        N = ECC*sinh(H) - H;
        N0 = N - sqrt(mu_val/((-SMA)^3))*dt;
        coord_set_out = [SMA,ECC,INC,RAAN,ARGP,wrapTo2Pi(N0)];
    else
        error('Code cannot convert parabolic orbit sets')
    end
    
elseif strcmpi(coord_type,'orb el') % Convert from orb el to cartesian
    
    SMA = coord_set(1);
    ECC = coord_set(2);
    INC = coord_set(3);
    RAAN = coord_set(4);
    ARGP = coord_set(5);
    
    alpha_check = 1/SMA;
    if alpha_check > 0 % Elliptical Case
        M0 = coord_set(6);
        M = M0 + sqrt(mu_val/(SMA^3))*dt;
        E = fzero(@(x) M - x + ECC*sin(x),M);
        TA = 2*atan(sqrt((1+ECC)/(1-ECC))*tan(E/2));
    elseif alpha_check < 0 % Hyperbolic case
        N0 = coord_set(6);
        N = N0 + sqrt(mu_val/((-SMA)^3))*dt;
        H = fzero(@(x) N + x - ECC*sinh(x),N);
        TA = 2*atan(sqrt((ECC + 1)/(ECC - 1))*tanh(H/2));
    else
        error('Code cannot convert parabolic orbit sets')
    end
    
    theta = ARGP + TA;
    r_mag = SMA*(1-ECC^2)/(1 + ECC*cos(TA));
    h_mag = sqrt(mu_val*SMA*(1-ECC^2));
    
    r_vec = r_mag*[cos(RAAN)*cos(theta) - sin(RAAN)*sin(theta)*cos(INC);...
                   sin(RAAN)*cos(theta) + cos(RAAN)*sin(theta)*cos(INC);...
                   sin(theta)*sin(INC)];
    
    v_vec = -mu_val/h_mag*...
        [cos(RAAN)*(sin(theta) + ECC*sin(ARGP)) + sin(RAAN)*(cos(theta) + ECC*cos(ARGP))*cos(INC);...
         sin(RAAN)*(sin(theta) + ECC*sin(ARGP)) - cos(RAAN)*(cos(theta) + ECC*cos(ARGP))*cos(INC);...
         -(cos(theta) + ECC*cos(ARGP))*sin(INC)];

     coord_set_out = [r_vec',v_vec'];
     
else
    error('coord_set argument must be either ''cart'' or ''orb el''')
end

end