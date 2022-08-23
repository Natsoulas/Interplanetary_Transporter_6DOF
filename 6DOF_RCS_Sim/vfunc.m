function v = vfunc(A,B,E,e,n)
    v = (n/(1 - e*cos(E)))*(-A*sin(E) + B*cos(E));
end