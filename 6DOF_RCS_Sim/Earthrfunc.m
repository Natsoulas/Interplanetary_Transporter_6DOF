function r = Earthrfunc(t,IC)
%Please replace Earth with desired body in function title,
%and keplerian elements with desired body's in the code below in order to include another body
%don't forget to update dynamics with gravity from the newly added body too
%should be orbit_attitude_numint.m or something similar for dynamics
    a = (1.495468303960526E+08)*1000; %km 2 meters
    e = 1.699803590732099E-02; % 
    I = 4.787963974336582E-03*pi/180;
    w = 2.700238217175390E+02*pi/180;
    O = 1.917131048127610E+02*pi/180;
    t_p = 2465791.336754327174; %JD5
    [A,B]= ABfunc(a,e,I,w,O);
    M = Mfunc(sqrt(IC.mu/a^3)*86400,t,t_p);
    E = invertKepler(M,e,M);
    r = rfunc(A,B,E,e);
end
