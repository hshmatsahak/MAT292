%  Details: Consider the second-order ODE
%
%  y'' + p(t) y' + q(t) y = g(t)
% 
% Write a second-order ODE solver using the method described above.
% 
% This m-file should be a function which accepts as variables 
% (t0,tN,y0,y1,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0, y1 are the initial conditions of
% the ODE, and h is the stepsize.  You may also want to pass the functions 
% into the ODE the way |ode45| does (check MATLAB lab 2). Name the function
% DE2_<UTORid>.m.
function second_order_solver = DE2_sahakhsh(p, q, g, t0, tN, y0, y1, h)
    num_steps = (tN-t0)/h;
    y = zeros(1, num_steps+1);
    t = t0:h:tN;
    
    y(1) = y0;
    y(2) = y(1) + y1*h;
    
    for i = 2:num_steps
        yp = (y(i) - y(i-1))/h;
        ypp = g(t(i)) -q(t(i))*y(i) - p(t(i))*yp;
        y(i+1) = h*h*ypp+2*y(i)-y(i-1);
    end
    
    second_order_solver = struct('t', t, 'y', y);
    
        
        
        
        