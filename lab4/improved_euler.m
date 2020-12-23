% Details: This m-file should be a function which accepts as variables 
% (t0,tN,y0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0 is the initial condition of the
% ODE, and h is the stepsize.  You may also want to pass the function into
% the ODE the way |ode45| does (check lab 2).

function improved_euler = improved_euler(f, t0, tN, y0, h)
num_steps = floor((tN - t0)/h);
y = zeros(1, num_steps+1);
y(1) = y0;
for i = 1:num_steps
    init_time = t0 + h*(i-1);
    init_y = y(i);
    end_time = t0 + h*(i);
    slope1 = f(init_time, init_y);
    end_y = init_y+h*slope1;
    slope2 = f(end_time, end_y);
    avg_slope = (slope1+slope2)/2;
    increment = avg_slope*h;
    y(i+1) = y(i) + increment;
end
x = linspace(t0, tN, num_steps+1);
improved_euler = struct('x', x, 'y', y);