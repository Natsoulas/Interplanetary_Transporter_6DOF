%RCS Thruster Placement Tool % more efficiet with its wrapper script
% written by Nicholas Natsoulas in 2022 under mentorship of Brent Faller @
% NASA Glenn Research Center
function Thruster_PlacementData = thrusterplacement(Thruster_PlacementData,noRings, config, heightarray, MOI, R, COM, defaultForce,move,dup,editForce,customthruster,visualize)
% Basic inputs - noRings: scalar value of number of rings to be placed
%       - config: which thruster configuration to use (provides number of
%       thrusters and their azimuthal locations with respect to Xb and
%       pointing diretions in cartesian body frame
%       - heightarray: what are the desired heights ex: [3, 2, 5]
%       - MOI: moment of inertia matrix
%       - R: Radius of Ring (bounds spacecraft bus)
%       - COM: Center of mass wrt to ring geometry (cartesian [x,y,z])
%       - defaultForce: Force for all thrusters (subject to editing later)
%
% Editing feature inputs organized by feature:
%      - move (move ring): needs parameters ring number - ringno and height - h
%        input architecture:  "ringid h" where string input is split at the
%        space, and h is converted into a double
%      - dup (duplicate ring): "ringid h" where string input is split at
%        the space, and h is converted into a double
%      - edit thruster force: needs parameters id, thruster number, force
%      - input architecture: "id thrusterno force" where string input is
%        split at the space, and force is converted into a double
%      - define custom thruster: needs parameters angle, height, Xpoint,
%        Ypoint, Zpoint, and Force
%      - input architecture: "angle h Xpoint Ypoint Zpoint Force" where all
%        parameters are split and then converted from string to a double
%
% Specific naming conventions for ids (important for parsing and optional
%      - inputs: rings ('ring1' 'ring2' 'customthruster1' customthruster2')
%      so that upon parsing string containing ring and or custom thruster
%      the tool knows to start a new input rather than concatenating it to
%      the previous input leading to an error. I.E. input: 'id angle
%      force', ' id angle force' as two inputs rather than 'id angle force
%      id angle force' which would create an error for whatever feature is
%      to be used.

%sample config input
%[thrustno,angle,Xpoint,Ypoint,Zpoint;thrustno2,angle,Xpoint,Ypoint,Zpoint;etc]

%if the Thruster_Placement Data is given a non-empty matrix then
%initialize:

if isempty(Thruster_PlacementData) == true

%initialize arrays subject to growth
Thruster_PlacementData = [];
rings = [];
customThrusters = [];
%%
%Main Tool actions/programming below
%ring initialization
for ringno = 1:noRings
rings = [rings; makeRing(ringno,config,heightarray(ringno),defaultForce)];
end
Thruster_PlacementData = [rings;customThrusters];
end
%optional editing programming below
if lower(move) ~= "none"
    %then parse with split function
    parameter_list = split(move);
    ringid = parameter_list(1);
    h = str2double(parameter_list(1));
    rings = moveRing(rings,ringid,h);
    Thruster_PlacementData = [rings;customThrusters];
end
if lower(dup) ~= "none"
    %then parse with split function
    param_list = split(dup);
    ringid = param_list(1);
    h = str2double(param_list(2));
    rings = dupRing(rings,ringid,h,noRings);
    Thruster_PlacementData = [rings;customThrusters];
end
if lower(editForce) ~= "none"
    %then parse with split function
    param_list = split(editForce);
    id = param_list(1);
    thrusterno = param_list(2);
    newforce = param_list(3);
    Thruster_PlacementData = editThrusterForce(Thruster_PlacementData,id,thrusterno, newforce);
end
if lower(customthruster) ~= "none"
    %then parse with spllt function and also check duplicates
    param_list = split(customthruster);
    id = param_list(1);
    angle = param_list(2);
    h = param_list(3);
    Xpoint = param_list(4);
    Ypoint = param_list(5);
    Zpoint = param_list(6);
    Force = param_list(7);
    customThruster = makecustomThruster(id,angle,h,Xpoint,Ypoint,Zpoint,Force);
    Thruster_PlacementData = [Thruster_PlacementData; customThruster];
end

%%
%visualization programming below
if visualize == true
    %placement visualization below
    figure
    hold on
    [A,B,C] = cylinder(R);
    hache = max(heightarray);
    C = C*hache;
    D = [];
    C = C -2;
    surf(A,B,C)
    for k = 1:1:size(Thruster_PlacementData,1)
        angle = str2double(Thruster_PlacementData(k,3));
        [X,Y] = azim2cartbody(angle,COM,R);
        h = str2double(Thruster_PlacementData(k,4));
        ZwrtCOM = h - COM(3);
        Z = ZwrtCOM;
        U = str2double(Thruster_PlacementData(k,5));
        V = str2double(Thruster_PlacementData(k,6));
        W = str2double(Thruster_PlacementData(k,7));
        force = str2double(Thruster_PlacementData(k,8));
        quiver3(X,Y,Z,U,V,W,force,'or')
        %textscatter3(X,Y,Z, Thruster_PlacementData(k,1))
    end
    %torque and angular momentum visualizations below
    %compute total torque vectors for all unique sets of n-thrusters
    dimlong = size(Thruster_PlacementData,1);
    thrusterchoices = 1:dimlong;
    combos = nchoosek(thrusterchoices,4);%unique combo ids %FIX THIS!!!
    ComboTorqs = [];
    for k = 1:1:size(combos,1)
        T = [0,0,0];
        for i = 1:1:size(combos,2)
            angle = str2double(Thruster_PlacementData(combos(k,i),3));
            [X,Y] = azim2cartbody(angle,COM,R);
            h = str2double(Thruster_PlacementData(combos(k,i),4));
            ZwrtCOM = h - COM(3);
            Z = ZwrtCOM;
            r = [X,Y,Z];
            U = str2double(Thruster_PlacementData(combos(k,i),5));
            V = str2double(Thruster_PlacementData(combos(k,i),6));
            W = str2double(Thruster_PlacementData(combos(k,i),7));
            force = str2double(Thruster_PlacementData(combos(k,i),8));
            F = force*[U,V,W];
            Tconcat = cross(r,F);
            T = T + Tconcat;
        end
        ComboTorqs = [ComboTorqs;T];
    end
    %Compute total angular velocities for all unique sets of n-thrusters
    % [angularvelocityfrom10milliseconds, angularvelocityfrom1sec]
    %angularvelocities =[];
    %for j = 1:1:size(ComboTorqs,1)
    %    angularacceleration = MOI\transpose(ComboTorqs(j,:));
    %    angularvelocity10mill = angularacceleration*0.01;
    %    angularvelocity1sec = angularacceleration*1;
    %    angularvelocities = [angularvelocities;transpose(angularvelocity10mill),transpose(angularvelocity1sec)];
    %end
    %convert cartesian quantities to spherical coordinates
    [aztrq,eltrq,rhotrq] = cart2sph(ComboTorqs(:,1),ComboTorqs(:,2),ComboTorqs(:,3));
    %[azavel10,elavel10,rhoavel10] = cart2sph(angularvelocities(:,1),angularvelocities(:,2),angularvelocities(:,3));
    %[azavel1,elavel1,rhoavel1] = cart2sph(angularvelocities(:,4),angularvelocities(:,5),angularvelocities(:,6));
    aztrq = rad2deg(wrapTo2Pi(aztrq)); eltrq = rad2deg(eltrq);
    %azavel10 = rad2deg(wrapTo2Pi(azavel10)); elavel10 = rad2deg(elavel10);
    %azavel1 = rad2deg(wrapTo2Pi(azavel1)); elavel1 = rad2deg(elavel1);
    %Locate all unique instances of [AZ EL] using [C,IA,IC] = unique([AZ,El],'rows')
    %[C,ia,ic] = unique([aztrq,eltrq],'rows');
    Unsorted = [aztrq,eltrq,rhotrq];
    Sorted = sortrows([aztrq,eltrq,rhotrq]);
    uniqueAzEl_MaxTorque = [];
    counter = 0;
    Sorted = sortrows([aztrq,eltrq,rhotrq]);
    uniqueAzEl_MaxTorque = [];
    counter = 0;
    for q = 1:1:size(Sorted,1)
        counter = counter + 1;
        %first condition must be that if it is not the same as the index
        %before and after it the row is therefore unique and appended for
        %later use
        if counter==1 && Sorted(q,1)~=Sorted(q+1,1) && Sorted(q,2) ~= Sorted(q+1,2)
            uniqueAzEl_MaxTorque = [uniqueAzEl_MaxTorque; Sorted(q,:)];
        elseif counter == size(Sorted,1)
            uniqueAzEl_MaxTorque = [uniqueAzEl_MaxTorque; Sorted(q,:)];
        %elseif Sorted(q,1)~=Sorted(q+1,1) && Sorted(q,2)~=Sorted(q+1,2) && Sorted(q,1)~=Sorted(q-1,1) && Sorted(q,2) ~= Sorted(q-1,2)
        %    uniqueAzEl_MaxTorque = [uniqueAzEl_MaxTorque; Sorted(q,:)];
        elseif (Sorted(q,1)~=Sorted(q+1,1) || Sorted(q,2)~=Sorted(q+1,2))%% && (Sorted(q,1)==Sorted(q-1,1)&&Sorted(q,2)==Sorted(q-1,2))
            uniqueAzEl_MaxTorque = [uniqueAzEl_MaxTorque; Sorted(q,:)];
        end
    end
    pointsize = 100;
    figure
    scatter(uniqueAzEl_MaxTorque(:,1),uniqueAzEl_MaxTorque(:,2),pointsize,uniqueAzEl_MaxTorque(:,3),'filled')
    xlim([0 360])
    ylim([-90 90])
    xlabel('Azimuth (Degrees)','FontSize',14)
    ylabel('Elevation (Degrees)','FontSize',14)
    c = colorbar;
    c.Label.String = 'Torque magnitude (N*m)';
    c.Label.FontSize = 14;
    title('\fontsize{18}Torque Scatterplot by Az & El')
    
    %
    angularvelocities = [];
    for j = 1:1:size(uniqueAzEl_MaxTorque,1)
    [x_o,y_o,z_o] = sph2cart(deg2rad(uniqueAzEl_MaxTorque(j,1)),deg2rad(uniqueAzEl_MaxTorque(j,2)),uniqueAzEl_MaxTorque(j,3));
    thing = [x_o;y_o;z_o];
    angularacceleration = MOI\thing;
    angularvelocity10mill = angularacceleration*0.01;
    angularvelocity1sec = angularacceleration*1;
    [a1,b1,c1] = cart2sph(angularvelocity10mill(1,1),angularvelocity10mill(2,1),angularvelocity10mill(3,1));
    angularvelocity10millazel = [a1,b1,c1];
    [a2,b2,c2] = cart2sph(angularvelocity1sec(1,1),angularvelocity1sec(2,1),angularvelocity1sec(3,1));
    angularvelocity1secazel = [a2,b2,c2];
    angularvelocities = [angularvelocities;angularvelocity10millazel,angularvelocity1secazel];
    end
    pointsize = 100;
    %10 sec burn scatterplot
    figure
    scatter(uniqueAzEl_MaxTorque(:,1),uniqueAzEl_MaxTorque(:,2),pointsize,angularvelocities(:,3),'filled')
    xlim([0 360])
    ylim([-90 90])
    xlabel('Azimuth (Degrees)','FontSize',14)
    ylabel('Elevation (Degrees)','FontSize',14)
    c = colorbar;
    c.Label.String = 'Angular Velocity magnitude (Degrees/s)';
    c.Label.FontSize = 14;
    title('\fontsize{18}Angular Velocity Scatterplot by Az & El (10 millisecond burn)')
    %1 sec burn scatterplot
    figure
    scatter(uniqueAzEl_MaxTorque(:,1),uniqueAzEl_MaxTorque(:,2),pointsize,angularvelocities(:,6),'filled')
    xlim([0 360])
    ylim([-90 90])
    xlabel('Azimuth (Degrees)','FontSize',14)
    ylabel('Elevation (Degrees)','FontSize',14)
    c = colorbar;
    c.Label.String = 'Angular Velocity magnitude (Degrees/s)';
    c.Label.FontSize = 14;
    title('\fontsize{18}Angular Velocity Scatterplot by Az & El (1 second burn)')
end
%%
%Define Subfunctions below
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
    function Ring = makeRing(ringno,config,h,defaultForce)
        Ring = [];
        for thrustno = 1:length(config(:,1))
            azimuthangle = config(thrustno,2);
            Xpoint = config(thrustno,3);
            Ypoint = config(thrustno,4);
            Zpoint = config(thrustno,5);
            ringid = strcat("ring ",num2str(ringno));
            Ring = [Ring;ringid,thrustno,azimuthangle,h,Xpoint,Ypoint,Zpoint,defaultForce];
        end
    end
    function rings = moveRing(rings,ringid,h)
        %moves a ring vertically along bus height
        %select desired ring by passing parameter for ring number (ringno)
        %pass second parameter for desired height to move to
        DesiredRing = [];
        desringidxs = [];
        for counter = 1:length(rings(:,1))
            if rings(counter,1) == ringid
                DesiredRing = [DesiredRing;rings(counter,1)];
                desringidxs = [desringidxs, counter];
            end
        end
        %update the heights of the rings within DesiredRing list
        DesiredRing(:,4) = h;
        %DesiredRing = [DesiredRing(:,1:3), h*ones(length(desringidxs),1),DesiredRing(:,5:end)];
        %switchout Old ring data with the updated DesiredRing
        if desringidxs(1) ~= 1 && desringidxs(end) < size(rings,1) %execute if desired ring is not the first ring, and is not the last ring
            rings = [rings(1:desringidxs(1)-1,:);DesiredRing;rings(1+desringidxs(end):end,:)];
        elseif desringidxs(1) == 1 && desringidxs(end) < size(rings,1) %execute if desired ring is the first ring but not the last ring
            rings = [DesiredRing;rings(1+desringidxs(end):end,:)];
        elseif desringidxs(1) ~= 1 && desringidxs(end) == size(rings,1)%execute if desired ring is not the first ring but is the last ring
            rings = [rings(1:desringidxs(1)-1,:);DesiredRing];
        else %execute if desired ring is the first ring and is the last ring
            rings = DesiredRing;
        end
    end
    function rings = dupRing(rings,ringid,h,noRings)
        %duplicates preexisting ring to a different inputted vertical
        %height (h)
        %other inputs - rings (info on all rings), ringno (ring number to be
        %duplicated)
        DesiredRing = [];
        desringidxs = [];
        %method for assigning desired heights to the desired ring data
        for counter = 1:length(rings(:,1))
            if rings(counter,1)==ringid
                DesiredRing = [DesiredRing;rings(counter,1)];
                desringidxs = [desringidxs, counter];
            end
        end
        DesiredRing(:,4) = h;
        %append DesiredRing as a new ring
            %first DesiredRing's assigned ring number (identifier) needs to be
            %reassigned to a new Ring number that isn't already in use
            %ring number that doesn't exist yet is preassinged number of rings
            %(noRings) plus 1, or (noRings + 1)
        DesiredRing(:,1) = noRings+1;
        rings = [rings;DesiredRing];
    end
    function Thruster_PlacementData = editThrusterForce(Thruster_PlacementData,ringid,thrustno, newForce)
        for rowidx = 1:size(Thruster_PlacementData,1)
            if Thruster_PlacementData(rowidx,1) == ringid && Thruster_PlacementData(rowidx,2) == thrustno
                Thruster_PlacementData(rowidx,8) = newForce;
            end
        end
    end
    function customThruster = makecustomThruster(customthrusterid,angle,h,Xpoint,Ypoint,Zpoint,Force)
        if nargin < 7
            Force = defaultForce;
        end
        %define custom thruster
        customThruster = [customthrusterid,1,angle,h,Xpoint,Ypoint,Zpoint,Force];
    end
end