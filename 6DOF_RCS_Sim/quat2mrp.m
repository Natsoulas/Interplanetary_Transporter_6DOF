function mrp = quat2mrp(quat)
%put shadow thingy in here
%output shadow set when dot product of itself exceeds 1
s1 = quat(2)/(1+quat(1));
s2 = quat(3)/(1+quat(1));
s3 = quat(4)/(1+quat(1));
mrp = [s1;s2;s3];
if dot(mrp,mrp) > 1
    s1 = -quat(2)/(1-quat(1));
    s2 = -quat(3)/(1-quat(1));
    s3 = -quat(4)/(1-quat(1));
    mrp = [s1;s2;s3];
end
end