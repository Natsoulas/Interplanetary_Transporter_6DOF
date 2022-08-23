# 6-DOF-Simulation---Orbit-Attitude-Tracking
This is a 6-DOF Spacecraft Orbit and Attitude Simulation that utilizes two controller (1 for orbit &amp; 1 for attitude), one simplex optimization method implementation, a PWPF modulator, and ode45 (RK 4-5) numerical integration to simulate a closed-loop control system for the given desired attitude, orbit, and initial conditions of just about any particular Space Mission the user chooses. The only attitude actuators are RCS (Reaction Thrusters).


The main file for user interface is the Control_Loop_ode45_Orbit.m file which takes all the spacecraft mass properties, thruster configuration, the initial attitude and orbit conditions/state and the time it takes place. The main while loop is also where the convergence criteria for the sim are stated. For now I have the sim set for orbital velocity vector attitude tracking, which lines up the body frame (long)z-axis of a monolithic spacecraft with the inertial (Heliocentric J-2000) orbital velocity direction during perigee passage (true anomaly -30 to 30 degrees).

However, anyone can set up this simulation with any example scenario desired. Just be careful to set up IC's right, tune the controllers, optimizer, and signal modulator, and add any celestial bodies of interest by duplicating functions of either earth or mars and replacing their keplerian elements with the desired body's keplerian elements.

Disclaimer: I did not develop/derive these lyapunov stable nonlinear control laws or this implementation of Simplex optimization method, some nice folks from CU Boulder and Matlab file exchange did :).
