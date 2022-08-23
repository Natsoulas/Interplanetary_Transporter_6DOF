function burn = timedburn(t_start,t_end,Force_B,t)
    if round(t) >= t_start && round(t) <= t_end
        burn = Force_B;
    else
        burn = [0;0;0];
    end
 end