%% Funwork 1
% Evan Greene 
% 
% 2020-01-30

%% Derivation of the equations of motion
% 
% The equations modeling the motion of a double inverted pendulum on a cart
% can be derived from the Lagrangian equations of motion. The equations of
% motion are given by
% 
% $$
% \frac{\partial d}{\partial dt} \left( \frac{\partial L}{\partial \dot{q}}
% \right) - \frac{\partial L}{\partial q} = f
% $$
% 
% Where $L$ is the Lagrangian fucntion $L = K - U$, the difference between
% kinetic and potential energy, $f$ is the non-conservative force acting
% on the system, and $q$ is the generalized coordinate used to describe the
% system.
% 
% For the case of the double-inverted pendulum, the state of the system is
% described by three generalized coordinates, $x$, $\theta_1$, and
% $\theta_2$. Therefore, there are three equations of motion.
% 
% To find find the equations of motion, we must first find the kinetic and
% potential energy of the system. The kinetic energy $K = K_0 + K_1 + K_2$,
% where $K_0$ is the kinetic energy of the cart
% and pendulum, $K_1$ is the kinetic energy of mass 1,
% and $K_2$ is the kinetic energy of mass 2. The value of
% $K_0$ is
% 
% $$ K_0 = \frac{1}{2}M \dot{x}^2 $$
% 
% While the value of $K_1$ is
% 
% $$
% \begin{array}{rcl}
% K_1 & = & \frac{1}{2}m_1 v_1^2 \\
% & = & \frac{1}{2} m_1 \left[ (\dot x + l_1 \dot \theta_1\cos\theta_1)^2
% + (l_1 \dot \theta_1 \sin \theta)^2 \right] \\
% & = & \frac{1}{2} m_1 \dot x^2 +
% \frac{1}{2} m_1 l_1^2 \dot \theta_1^2 +
% m_1 l_1 \dot x \dot \theta_1 \cos \theta_1
% \end{array}
% $$
% 
% and the value of $k_2$ is
% 
% $$
% \begin{array}{rcl}
% K_2 & = & \frac{1}{2}m_2 v_2^2 \\
% & = & \frac{1}{2} m_2 \left[ (\dot x + l_1 \dot \theta_1 \cos \theta_1 +
% l_2 \dot{\theta}_2 \cos \theta_2)^2
% + (-l_1 \dot \theta_1 \sin \theta_1 -
% l_2 \dot{\theta}_2 \sin \theta_2)^2 \right] \\
% & = & \frac{1}{2} m_2 \dot x^2 + \frac{1}{2} l_1^2 \dot \theta_1^2 +
% \frac{1}{2} m_2 l_2^2 \dot \theta_2^2 +
% m_2 l_1 \dot x \dot \theta_1 \cos \theta_1 +
% m_2 l_2 \dot x \dot \theta_2 \cos \theta_2 \\ &  & \quad + \,
% m_2 l_1 l_2 \dot \theta_1 \dot \theta_2 \cos (\theta_2 - \theta_1)
% \end{array}
% $$
% 
% Meanwhile, the only form of potential energy present in the system is
% gravitational potential energy, which is associated with the height of
% each of the masses.
% 
% $$ U = m_1 g l_1 \cos \theta_1 +
% m_2 g \left( l_1 \cos \theta_1 + l_2 \cos \theta_2 \right) $$
% 
% The Lagrangian function is therefore
% 
% $$
% \begin{array}{rcl}
% L & = & \frac{1}{2} \left(M + m_1 + m_2\right) \dot{x}^2
% + \frac{1}{2} (m_1 + m_2) l_1^2 \dot{\theta}_1^2
% + \frac{1}{2} m_2 l_2^2 \dot{\theta}_2^2
% \\&& \quad
% + (m_1 + m_2) l_1 \dot x \dot \theta_1 \cos \theta_1
% + m_2 l_2 \dot x \dot \theta_2 \cos \theta_2
% + m_2 l_1 l_2 \dot \theta_1 \dot \theta_2 \cos (\theta_2 - \theta_1)
% \\&& \quad
% - (m_1 + m_2) g l_1 \cos \theta_1 - m_2 g l_2 \cos \theta_2
% \end{array}
% $$
% 
% We can now find the equations of motion from the Lagrangian. The first
% equation is
% 
% $$
% \frac{\partial d}{\partial dt} \left( \frac{\partial L}{\partial \dot{x}}
% \right) - \frac{\partial L}{\partial x} = u
% $$
% 
% We evaluate this function by finding
% 
% $$
% \frac{\partial L}{\partial \dot{x}} = (M + m_1 + m_2) \dot x
% + (m_1 + m_2) l_1 \dot \theta_1 \cos \theta_1
% + m_2l_2 \dot \theta_2 \cos \theta_2
% $$
% 
% $$
% \begin{array}{rcl}
% \frac{d}{dt}\frac{\partial L}{\partial \dot{x}} & = &(M + m_1 + m_2)\ddot{x}
% + (m_1 + m_2) l_1 \ddot \theta_1 \cos \theta_1
% + m_2 l_2 \ddot \theta_2 \cos \theta_2
% \\&& \quad
% - (m_1 + m_2) l_1 \dot \theta_1^2 \sin \theta_1
% + m_2 l_2 \dot \theta_2^2 \sin \theta_2
% \end{array}
% $$
% 
% $$
% \frac{\partial L}{\partial x} = 0
% $$
% 
% Therefore, the first equation of motion is
% 
% $$
% \begin{array}{rcl}
% u & = &
% (M + m_1 + m_2)\ddot{x}
% + (m_1 + m_2) l_1 \ddot \theta_1 \cos \theta_1
% + m_2 l_2 \ddot \theta_2 \cos \theta_2
% \\&& \quad
% - (m_1 + m_2) l_1 \dot \theta_1^2 \sin \theta_1
% + m_2 l_2 \dot \theta_2^2 \sin \theta_2
% \end{array}
% $$
% 
% We can find the second equation of motion by finding
% 
% $$
% \begin{array}{rcl}
% \frac{\partial L}{\partial \dot{\theta}_1} & = &
% (M + m_1 + m_2)\dot{x}
% + (m_1 + m_2) l_1 \dot{\theta}_1
% + (m_1 + m_2) l_1 \dot \theta_1 \cos \theta_1
% \\&& \quad
% + m_2 l_2 \dot \theta_2  \cos (\theta_2 - \theta_1)
% \end{array}
% $$
% 
% $$
% \begin{array}{rcl}
% \frac{d}{dt} \frac{\partial L}{\partial \dot{\theta}_1} & = &
% (m_1 + m_2) l_1 \ddot x \cos \theta_1
% + (m_1 + m_2) l_1^2 \ddot \theta_1
% + m_2 l_1 l_2 \ddot \theta_2 \cos(\theta_2 - \theta_1 )
% \\&& \quad
% - m_2 l_1 l_2 \dot \theta_2^2 \sin(\theta_2 - \theta_1)
% - (m_1 + m_2) \dot x \dot \theta_1 \sin \theta_1
% \\&& \quad
% + m_2 l_1 l_2 \dot \theta_1 \dot \theta_2 \sin(\theta_2 - \theta_1)
% \end{array}
% $$
% 
% $$
% \begin{array}{rcl}
% \frac{\partial L}{\partial \theta_1} & = &
% -(m_1 + m_2) \dot x \dot \theta_1 \sin \theta_1
% + m_2 l_1 l_2 \dot \theta_1 \dot \theta_2 \sin (\theta_2 - \theta_1)
% \\&& \quad
% - (m_1 + m_2) g l_1 \sin \theta_1
% \end{array}
% $$
% 
% The second equation of motion is therefore
% 
% $$
% \begin{array}{rcl}
% 0 &=&
% (m_1 + m_2) l_1 \ddot x \cos \theta_1
% + (m_1 + m_2) l_1^2 \ddot \theta_1
% + m_2 l_1 l_2 \ddot \theta_2 \cos(\theta_2 - \theta_1 )
% \\&& \quad
% - m_2 l_1 l_2 \dot \theta_2^2 \sin(\theta_2 - \theta_1)
% + (m_1 + m_2) g l_1 \sin \theta_1
% \end{array}
% $$
% 
% We can find the third equation of motion from
% 
% $$
% \begin{array}{rcl}
% \frac{\partial L}{\partial \dot{\theta}_2} & = &
% m_2 l_2^2 \dot \theta_2
% + m_2 l_2 \dot x \cos \theta_2
% + m_2 l_1 l_2 \dot \theta_2 \cos (\theta_2 - \theta_1)
% \end{array}
% $$
% 
% $$
% \begin{array}{rcl}
% \frac{d}{dt} \frac{\partial L}{\partial \dot{\theta}_2} & = &
% m_2 l_2 \ddot x \cos \theta_2
% + m_2 l_1 l_2 \ddot \theta_1 \cos(\theta_2 - \theta_1)
% + m_2 l_2 \ddot \theta_2
% \\&& \quad
% + m_2 l_1 l_2 \dot \theta_1^2 \sin(\theta_2 - \theta_1)
% - m_2 l_2 \dot x \dot \theta_2 \sin \theta_2
% \\&& \quad
% - m_2 l_1 l_2 \dot \theta_1 \dot \theta_2 \sin (\theta_2 - \theta_1)
% \end{array}
% $$
% 
% $$
% \begin{array}{rcl}
% \frac{\partial L}{\partial \theta_1} & = &
% - m_2 l_2 \dot x \dot \theta_2 \sin \theta_2
% - m_2 l_1 l_2 \dot \theta_1 \dot \theta_2 \sin (\theta_2 - \theta_1)
% + m_2 g l_2 \sin \theta_2
% \end{array}
% $$
% 
% The third equaton of motion is therefore
% 
% $$
% \begin{array}{rcl}
% 0 & = &
% m_2 l_2 \ddot x \cos \theta_2
% + m_2 l_1 l_2 \ddot \theta_1 \cos(\theta_2 - \theta_1)
% + m_2 l_2 \ddot \theta_2
% \\&& \quad
% + m_2 l_1 l_2 \dot \theta_1^2 \sin(\theta_2 - \theta_1)
% - m_2 g l_2 \sin \theta_2
% \end{array}
% $$
% 
%% Derivation of the state-space model
% From (Bogdanov 2004) we can rearrange these equations into the form
% 
% $$
% \mathbf D (\mathbf q ) \ddot q + \mathbf C (\mathbf q, \mathbf{\dot q})
% + \mathbf g(\mathbf q) = \mathbf H u
% $$
% 
% where
% 
% $$
% \mathbf q = \left[\begin{array}{c}
% x \\ \theta_1 \\ \theta_2
% \end{array} \right]
% $$
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
% and
% 
% $$ \mathbf H = \left[ \begin{array}{c} 1 \\ 0 \\ 0 \end{array} \right] $$
% 
% We can transform this model into the state-space format by creating the state
% vector
% 
% $$
% \mathbf x =
% \left[ \begin{array}{c} q \\ \dot q \end{array} \right] =
% \left[\begin{array}{c} x \\ \theta_1 \\ \theta_2 \\
% \dot x \\ \dot \theta_1 \\ \dot \theta_2 \end{array}\right]
% $$
% 
% From our above equations, we can see that
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
% \mathbf 0 \\ -\mathbf D^{-1} \mathbf H
% \end{array} \right] u
% $$
%
% Thus, we have a state-space model in the form $ \dot x = f(x,u)$ 
% 
%% MATLAB model creation 
% To create our model, we start by establishing our constants
clear all; 
M = 1.5; % kg 
m1 = 0.5; % kg 
m2 = 0.75; % kg 
l1 = 0.5; % m
l2 = 0.75; % m 
g = 9.81; % m/s^2
m12 = m1 + m2; 
degrees = 180 / pi; % degrees / radian
% D, C, and G arrays
D = @(x) [M + m12, ...
        m12 * l1 * cos(x(2)), ...
        m2 * l2 * cos(x(3));
        m12* l1 * cos(x(2)), ...
        m12 * l1^2, ...
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
%% Linearization
% To design a linearized model of the plant, we create a symbolic function
% equivalent to the handle function established in the previous section. 
% This allows us to use the |jacobian| function to linearize. 
% The Taylor linearization of $\dot x = f(x, u)$ about a point 
% $(x_0, u_0)$, 
% is given by 
% 
% $$ 
% \mathbf{\dot x} \approx \mathbf f(\mathbf x_0, u_0)
% + \frac{\partial \mathbf f}{\partial \mathbf x} (\mathbf x - \mathbf x_0)
% + \frac{\partial \mathbf f}{\partial u} (u - u_0)
% $$
% 
% where $\frac{\partial \mathbf f}{\partial \mathbf x}$ and 
% $\frac{\partial \mathbf f}{\partial u}$ are the jacobian matrices of $f$.
% (Zak 2003 p.75) 
% Also, since $f(\mathbf x_0, u_0) = \mathbf 0$, we can remove this term 
% from the equation. 

% create symbolic variables so the Jacobian function will work.
syms x [6 1] 
syms u 

% Calculate the jacobians of f with respect to 
J_f_x = jacobian(f(x, u), x); 
J_f_u = jacobian(f(x, u), u);

% evaluate those jacobians at x = 0 and u = 0 to get our linearization. 
A_linear = double(subs(J_f_x, x, zeros(6, 1)))
b_linear = double(subs(J_f_u, x, zeros(6, 1)))

%% Controllability and observability
% 
% If we consider only the single-input continuous case, a system 
% 
% $\mathbf{\dot x} = \mathbf{A} \mathbf x + \mathbf b u$
% 
% is controllable if and only if the controllability matrix 
% 
% $$
% rank \left[ \begin{array}{cccc}
% \mathbf b & \mathbf{Ab} & \cdots & \mathbf A^{n-1} \mathbf b
% \end{array} \right] = n
% $$
% 
% where $n$ is the size of the square matrix $\mathbf A$. (Zak 2003, p.95) 
% 
% This matrix can be found with 
controllable = rank(ctrb(A_linear, b_linear)) == 6

%%
% Since the controllability matrix is full rank the double-inverted 
% pendulum on 
% a cart is controllable. 

%% 
% The observability matrix is similar in construction. A multi-
% output system 
% 
% $$ \mathbf{\dot x} = \mathbf{Ax} + \mathbf b u$$ 
% 
% $$ \mathbf y = \mathbf{Cx} + \mathbf{D}u $$
% 
% or equivalently the pair $(\mathbf A, \mathbf C)$ is observable iff the 
% observability matrix 
% 
% $$ 
% rank
% \left[ \begin{array}{c}
% \mathbf C  \\ \mathbf{CA} \\ \vdots \\ \mathbf{CA}^{n-1}
% \end{array} \right] = n
% $$ 
% 
% where $n$ is the size of $A$.  (Zak 2003, p110)
% 
% To determine the observability of our system we must first find the 
% relationship between state and output of our model. 
% 
% $$ \mathbf y = \mathbf C \mathbf x $$ 
% 
% The $\mathbf C$ here should not be confused with the $\mathbf C$ that 
% corresponds to the coriolis and centrifugal forces in our equations of 
% motion. 
% 
% In this case the output is simply the first three elements of the state
% vector, so 
C_linear = [eye(3), zeros(3)];

%% 
% We can then test for observability

observable = rank(obsv(A_linear, C_linear)) == 6
%%
% Since the observability matrix has full column rank, the linearized 
% double-inverted pendulum on cart is observable. 
%% Transfer Function

% The transfer function can then be computed with 
state_space_linear = ss(A_linear, b_linear, C_linear, zeros(3, 1)); 
transer_function_linear = tf(state_space_linear)


%% Animation
% close all open figures 
close all
figure(1)

% Set the paramters for Euler integration
tfinal = 10; % seconds of animation
dt = 0.001; % step size
time = linspace(0, tfinal, tfinal / dt);

% set initial conditions
x_init = [0 0.01 0.02 0 0 0]';
x_current = x_init;

% Set up arrays for logging data here 
x_log = zeros(length(x_init), length(time));

% Create graphical elements 
% Rotation matrices
R1 = @(x) [cos(x(2)), -sin(x(2)); sin(x(2)), cos(x(2))];
R2 = @(x) [cos(x(3)), -sin(x(3)); sin(x(3)), cos(x(3))];
 

% the location of the first mass from the base. 
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
floor = line('xdata', [-2, 2], ...
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
% the step size is very small, so only update the graphics like once every
% few frames. 
frameRate = 60;
stepsPerFrame = 1 / (dt * frameRate);

% flag for whether to record video.
record_video = 0;

% if the flag is true, create the movie
if (record_video)
    % create movie 
    v = VideoWriter('DIPC.avi');
    v.FrameRate = frameRate; 
    v.open()
end

% perform the euler integration. 
for index = 1:length(time) - 1
    % find input as a function of the state in the last time step
    u = 0;  % no controller for now. 
    
    % update plant model. 
    dx_current = f(x_current, u); % find dx/dt
    % Euler integration -- x[k] = x[k-1] + dx[k-1]/dt * dt
    x_current = x_current + dx_current * dt;
    
    % Perform logging here
    x_log(:, index) = x_current;
    
    % allows the fps of the animation to be different from the euler 
    % integration step size. 
    if mod(index, stepsPerFrame) < .999
        % update point1 and point2
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
        if (record_video)
            frame = getframe;
            writeVideo(v, frame);
        end
    end 
end 

if (record_video)
    close(v);
end

%% 
% some post-mortem analysis
figure(2)

% plot position vs time 
subplot(231)
plot(time, x_log(1, :))
xlabel('t'); ylabel('x(t)');

% plot theta_1 vs time 
subplot(232)
plot(time, x_log(2, :))
xlabel('t'); ylabel('$\theta_1 (t)$', 'interpreter', 'latex');

% plot(theta_2) vs time 
subplot(233)
plot(time, x_log(3, :))
xlabel('t'); ylabel('$\theta_2 (t)$', 'interpreter', 'latex');

% plot dx/dt vs time 
subplot(234)
plot(time, x_log(4, :))
xlabel('t'); ylabel('$\frac{dx}{dt}$', 'interpreter', 'latex');

% plot d\theta_1/dt vs time 
subplot(235) 
plot(time, x_log(5, :))
xlabel('t'); ylabel('$\frac{d\theta_1}{dt}$', 'interpreter', 'latex');

% plot d\theta_2/dt vs time 
subplot(236) 
plot(time, x_log(6, :))
xlabel('t'); ylabel('$\frac{d\theta_2}{dt}$', 'interpreter', 'latex');


%% Bibliography
% Bogdanov, Alexander. ?Optimal Control of a Double Inverted Pendulum on 
% a Cart. Technical Report CSE-04-006, December 2004. 
% 
% Zak, Stanislaw H. Systems and Control. Oxford University Press, 2003.


    
