%% Funwork 3
%
%%  Evan Greene 
% 
% 2020-03-05
clear;
close all; 

%% Creating a DIPC model
%
% We use the same DIPC model that was created in Funwork 1 and used in 
% Funwork 2. That model is given by the equation
% 
% $$
% \mathbf{\dot x} =
% \left[ \begin{array}{cc}
% \mathbf 0 & \mathbf I_3 \\
% \mathbf 0 & -\mathbf D^{-1} \mathbf C
% \end{array} \right] x
% +
% \left[ \begin{array}{c}
% \mathbf 0 \\ -\mathbf D^{-1} \mathbf g
% \end{array} \right]
% +
% \left[ \begin{array}{c}
% \mathbf 0 \\ \mathbf D^{-1} \mathbf H
% \end{array} \right] u
% $$
%
% 
% where 
% $$ x = \left[ \begin{array}{cccccc} 
% x & \theta_1 & \theta_2 & \dot x & \dot \theta_1 & \dot \theta_2 
% \end{array} \right] $$
% 
% $$
% \mathbf D = \left[ \begin{array}{ccc}
% M + m_1 + m_2 &
% (m_1 + m_2) l_1  \cos \theta_1 &
% m_2 l_2  \cos \theta_2 \\
% (m_1 + m_2) l_1 \cos \theta_1 &
% (m_1 + m_2) l_1^2 &
% m_2 l_1 l_2  \cos(\theta_2 - \theta_1 ) \\
% m_2 l_2  \cos \theta_2 &
% m_2 l_1 l_2  \cos(\theta_2 - \theta_1) &
% m_2 l_2^2
% \end{array} \right]
% $$
% 
% $$
% \mathbf C = \left[ \begin{array}{ccc}
% 0 &
% - (m_1 + m_2) l_1 \dot \theta_1 \sin \theta_1 &
% m_2 l_2 \dot \theta_2 \sin \theta_2 \\
% 0 &
% 0 &
% - m_2 l_1 l_2 \dot \theta_2 \sin(\theta_2 - \theta_1) \\
% 0 &
% - m_2 l_1 l_2 \dot \theta_1 \sin(\theta_2 - \theta_1) &
% 0
% \end{array} \right]
% $$
% 
% $$
% \mathbf g = \left[ \begin{array}{c}
% 0 \\
% (m_1 + m_2) g l_1 \sin \theta_1 \\
% m_2 g l_2 \sin \theta_2
% \end{array}\right]
% $$
% 
% $$ \mathbf H = \left[ \begin{array}{c} 1 \\ 0 \\ 0 \end{array} \right] $$
% 
% Or, in MATLAB code

% Establish constants
M = 1.5; % kg 
m1 = 0.5; % kg 
m2 = 0.75; % kg 
l1 = 0.5; % m
l2 = 0.75; % m 
g = 9.81; % m/s^2
m12 = m1 + m2; 
degrees = 180 / pi; % degrees / radian

D = @(x) [M + m12, ...
        m12 * l1 * cos(x(2)), ...
        m2 * l2 * cos(x(3));
        m12* l1 * cos(x(2)), ...
        m12* l1^2, ...
        m2 * l1 * l2 * cos(x(3) - x(2));
        m2 * l2 * cos(x(3)), ...
        m2 * l1 * l2 * cos(x(3) - x(2)), ...
        m2 * l2^2 ];

C = @(x) [ 0, -m12 * l1 * x(5) * sin(x(2)), -m2 * l2 * x(6) * sin(x(3));
        0, 0,             -m2 * l1 * l2 * x(6) * sin(x(3) - x(2));
        0, m2 * l1 * l2 * x(5) * sin(x(3) - x(2)), 0];

G = @(x) [0,-m12 * g * l1 * sin(x(2)), -m2 * g * l2 * sin(x(3))]';

H = [1 0 0]';

f = @(x, u) [zeros(3),  eye(3); ...
            zeros(3),   -D(x) \ C(x)] * x ...
          + [zeros(3, 1); ...
            -D(x)\G(x)] ...
            + [zeros(3, 1); ...
              D(x)\H]*u;

%% Existance of a non-zero equilibrium for a one-input DIPC
% 
% For any time-invariant system such as our DIPC, the pair $(x_e,\, u_e)$ 
% is an equilibrium if $\dot x  = f(x_e,\, u_e) = 0$. 
% It is possible to show that for any non-zero state $x_e$, there is no 
% input $u_e$ that will bring the system into equilibrium. 
% 
% $$ x_e = \left[ \begin{array}{cccccc} 
% 0.1 & 60^\circ & 45^\circ & 0 & 0 & 0
% \end{array} \right] $$ 
% 
% we find
 
x_e = [0.1,  deg2rad(60),    deg2rad(45), 0,  0,  0]';
syms u 

disp("f(x_e, u) = ")
vpa(f(x_e, u), 4)
%% 
% We can see by inspection that the system of equations is overconstrained 
% and that no solution exists,
% but a more rigorous way of demonstrating this is 
% the |solve| function

disp("If the answer is an empty sym, there is no equilibrium")
solve(f(x_e, u) == zeros(1, 6)) 

%% Existence of a non-zero equilibrium for a two-input DIPC
%
% A two-input DIPC has the same problem as a first. In solving the equation
% $f(x_e, u) = 0$ for $u$, there are two unknowns and three equations, so 
% the system of equations is overconstrained and no solution exists. 

% modify the function f for a two input DIPC. 
% Replace H and update f to reflect the new H;
H = [1  0;
    0   1;
    0   0;];

f = @(x, u) [zeros(3),  eye(3); ...
            zeros(3),   -D(x) \ C(x)] * x ...
          + [zeros(3, 1); ...
            -D(x)\G(x)] ...
            + [zeros(3, 2); ...
              D(x)\H]*u;

% change the dimensionality of $u$
syms u [2 1]

%% 
% The equation that results is 
disp("f(x_e, u) = ")
vpa(f(x_e, u), 4)

%%
% and the solution is 
disp("If the answer is an empty sym, there is no equilibrium")
solution = solve(f(x_e, u) == zeros(1, 6)); 
solution.u1
solution.u2

%% Existence of a non-zero equilibrium for a three-input DIPC 
%
% From the previous two sections, it follows that the equations 
% $f(x_e, u = 0)$ will need three inputs for a solution to exist. 
% If we create a three-input system 
H = eye(3); 

f = @(x, u) [zeros(3),  eye(3); ...
            zeros(3),   -D(x) \ C(x)] * x ...
          + [zeros(3, 1); ...
            -D(x)\G(x)] ...
            + [zeros(3); ...
              D(x)\H]*u;

% change the dimensionality of u
syms u [3 1]

%% 
% The equation becomes
disp("f(x_e, u) = ")
vpa(f(x_e, u), 4)

%% 
% and the solution is 
solution = solve(f(x_e, u) == zeros(1, 6)); 
u_e = double([solution.u1; solution.u2; solution.u3])

%% Linearization of the three-input DIPC about the non-zero equilibrium
% The Taylor linearization of $\dot x = f(x, u)$ about a point 
% $(x_e, u_e)$, 
% is given by 
% 
% $$ 
% \mathbf{\dot x} \approx \mathbf f(\mathbf x_e, u_e)
% + \frac{\partial \mathbf f}{\partial \mathbf x} (\mathbf x - \mathbf x_e)
% + \frac{\partial \mathbf f}{\partial u} (u - u_e)
% $$
% 
% where $\frac{\partial \mathbf f}{\partial \mathbf x}$ and 
% $\frac{\partial \mathbf f}{\partial u}$ are the jacobian matrices of $f$.
% (Zak 2003 p.75) 
% Also, since $f(\mathbf x_0, u_0) = \mathbf 0$, we can remove this term 
% from the equation. 
%
% We can use MATLAB's symbolic variables and |jacobian| function to compute
% the Jacobian matrices. 

% create symbolic variables so the Jacobian function will work.
syms x [6 1] 
syms u [3 1]

% Calculate the jacobians of f with respect to 
J_f_x = jacobian(f(x, u), x); 
J_f_u = jacobian(f(x, u), u);

% evaluate those jacobians at x = x_e and u = u_e to get our linearization. 
A_linear = double(subs(subs(J_f_x, x, x_e), u, u_e))
b_linear = double(subs(subs(J_f_u, u, u_e), x, x_e))

clear x u

%% Controller Design Using Linear Matrix Inequalties

% TODO

%% Animation
% To easily animate the DIPC, we can create a function to save having to
% repeat work. 

% we can check whether this function works by calling animate with 
% the controller set to zero and the observer set to the actual state. 
controller = @(x) zeros(3, 1);
observer = @(xhat, y, u)  zeros(6, 1);
x_init = [0 0.01 0.02 0 0 0]';

animate(f, controller, observer, x_init)


function animate(f, controller, observer, x_init)
% inputs -- 
% f             -   the state-space function dx/dt  = f(x, u)
% controller    -   the function that relates the control input to the 
%                   estimated state, u = controller(~x)
% observer      -   the function that relates the input and output to 
%                   the estimated state d~x/dt = observer(~x, x, u)
% x_init        -   the initial state of the system.
% outputs -- 
% None. Plays animation
% close all open figures 
close all
figure(1)

% Set the paramters for Euler integration
tfinal = 5; % seconds of animation

% create the system of equations for the ODE solver to solve. 
odefun = @(t, augmented_state) ...
    [f(augmented_state(1:6, :), controller(augmented_state(7:12, :)));
    observer(augmented_state(7:end, :), augmented_state(1:3, :), ...
    controller(augmented_state(1:6, :)))];
% function d_augmented_state = odefun(t, augmented_state)
%     % inputs -- 
%     % t                 -   the time input
%     % augmented_state   -   the concatenated state variables x (actual
%     %                       state), ~x (estimated state)
%     % outputs -- 
%     % d_augmented_state -   the concatenated state variables dx (actual
%     %                       state) d~x
%     % unpack the augmented state
%     x = augmented_state(1:6);
%     x_est = augmented_state(7:end);
% 
%     % create the differential agumented state
%     dx = f(x, controller( x_est));
%     y = x(1:3);
%     dx_estimated = observer(x_est, y, controller(x_est));
% 
%     % pack the d_augmented_state
%     d_augmented_state = [dx; dx_estimated];
% end 
% create the augmented state
initial_augmented_state = [x_init; zeros(size(x_init))]; 
% solve the ODE 
[time, augmented_state] = ode45(odefun, [0, tfinal], ...
    initial_augmented_state);
% unpack the augmented state
time = time'; % make time a row vector
x_actual = augmented_state(:, 1:6)'; % make x_current a row vector
x_estimated = augmented_state(7:12)'; % make x_estimated a row vector

% animate the solution
% Create graphical elements 
% Define everything with respect to the first point
x_current = x_actual(:, 1);
% Rotation matrices
R1 = @(x) [cos(x(2)), -sin(x(2)); sin(x(2)), cos(x(2))];
R2 = @(x) [cos(x(3)), -sin(x(3)); sin(x(3)), cos(x(3))];

% the location of the first mass from the base. 
l1 = 0.5; l2 = 0.75; % m 
point1 = @(x) [x(1); 0] + R1(x) * [0; l1];
point1_current = point1(x_current);

% the location of the second mass 
point2 = @(x) point1(x) + R2(x) * [0;l2];
point2_current = point2(x_current);

% The size of the cart
cart_width = 1; cart_height = 0.25;
cart_position = [x_current(1) - 0.5*cart_width, -cart_height, ...
                cart_width, cart_height];
            
% a line for the floor 
ground = line('xdata', [-2, 2], ...
             'ydata', [-cart_height, -cart_height], ...
             'linewidth', 2, 'color', 'k');

% a rectangle for the cart
cart = rectangle('Position', cart_position, ... 
                    'EdgeColor', 'b', 'linewidth', 2);

% the hinge of the pendulum base
mass0 = line('xdata', double(x_current(1)), ...
             'ydata', 0, ...
             'linewidth', 3, 'color', 'r', 'marker', '*'); 

% line connecting the hinge and the first mass
bar1 = line('xdata', [x_current(1), point1_current(1)], ...
            'ydata', [0,            point1_current(2)], ...
            'linewidth', 2, 'color', 'b');
        
% the first mass object. 
mass1 = line('xdata', point1_current(1), ...
                'ydata', point1_current(2), ...
                'linewidth', 5, 'color', 'r', 'marker', '*');

% line connecting first and second masses
bar2 = line('xdata', [point1_current(1), point2_current(1)],...
            'ydata', [point1_current(2), point2_current(2)], ...
            'linewidth', 2, 'color', 'b');
        
% second mass 
mass2 = line('xdata', point2_current(1), ...
             'ydata', point2_current(2), ...
             'linewidth', 3, 'color', 'r', 'marker', '*');
            
% graph settings
axis([-2 2, -1.5, 1.5])
set(gca, 'dataaspectratio', [1 1 1])
axis on
grid on 
box on

% figure out the frame rate. 
% the step size is very small, so only update the graphics once every
% few frames. 
frameRate = 60;

index = 1;
dt = diff(time);
while index < length(time)
    % update point1 and point2
    x_current = x_actual(:, index);
    point1_current = point1(x_current);
    point2_current = point2(x_current);
    % set all the graphical elements.
    cart_position = [x_current(1) - 0.5*cart_width, -cart_height, ...
                    cart_width, cart_height];
    set(cart,   'Position', cart_position);
    set(mass0,  'xdata', x_current(1),              ...
                'ydata', 0);
    set(bar1,   'xdata', [x_current(1), point1_current(1)],   ...
                'ydata', [0,            point1_current(2)]);
    set(mass1,  'xdata', point1_current(1), ...
                'ydata', point1_current(2));
    set(bar2,   'xdata', [point1_current(1), point2_current(1)],...
                'ydata', [point1_current(2), point2_current(2)]);
    set(mass2,  'xdata', point2_current(1), ...
                'ydata', point2_current(2))  
    drawnow;
    % only draw at 60 fps, regardless of step size. 
    % this keeps the time steps consistent. 
    index = index + max([1, floor((1 / frameRate) / dt(index))]);
    
end 
% plot actual and estimated states
figure(2)
    for index = 1:6
        subplot(2, 3, index)
        plot(time, x_actual(index, :))
        hold on
        plot(time, x_estimated(index, :))
        hold off
        legend("x", "\hat{x}")
        xlabel("time (s)");
        ylabel(sprintf("x%d", index)); 
    end 
end 


