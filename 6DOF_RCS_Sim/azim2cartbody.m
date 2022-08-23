function [X,Y] = azim2cartbody(theta,COM,R)
        Xcom = COM(1);
        Ycom = COM(2);
        Xrel= R*cos(deg2rad(theta));
        Yrel = R*sin(deg2rad(theta));
        XwrtCOM = -Xcom + Xrel;
        YwrtCOM = -Ycom + Yrel;
        X = XwrtCOM;
        Y = YwrtCOM;
end