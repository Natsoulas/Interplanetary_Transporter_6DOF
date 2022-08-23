function r = Marsrfunc(t,IC)
    a = (2.279419507091813E+08)*1000;% km 2 meters
    e = 9.358077935297598E-02;
    I = deg2rad(1.846477690234122E+00); %degrees input
    w = deg2rad(2.867965399000748E+02);%
    O = deg2rad(4.943608223489434E+01);
    t_p = 2465934.922529480420;%in JD number figure out how to convert to seconds
    [A,B]= ABfunc(a,e,I,w,O);
    M = Mfunc(sqrt(IC.mu/a^3)*86400,t,t_p);
    E = invertKepler(M,e,M);
    r = rfunc(A,B,E,e);
end