%% Funwork 1
% Evan Greene 
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
% - (m_1 + m_2) g l_1 \cos \theta_1 + m_2 g l_2 \cos \theta_2
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
% + (m_1 + m_2) g l_1 \sin \theta_1
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
% - m_2 g l_2 \sin \theta_2
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
% + m_2 g l_2 \sin \theta_2
% \end{array}
% $$
% 
%% Derivation of the state-space model
% We can rearrange these equations into the form
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

% create a state vector 
syms x [6 1];
D = [M + m12, ...
        m12 * l1 * cos(x(2)), ...
        m2 * l2 * cos(x(3));
        m12* l1 * cos(x(2)), ...
        m12 * l1^2, ...
        m2 * l1 * l2 * cos(x(3) - x(2));
        m2 * l2 * cos(x(3)), ...
        m2 * l1 * l2 * cos(x(3) - x(2)), ...
        m2 * l2^2 ];
 
 C = [0, -m12 * l1 * x(5) * sin(x(2)), m2 * l2 * x(6) * sin(x(3));
     0, 0, -m2 * l1 * l2 * x(6) * sin(x(3) - x(2));
     0, -m2 * l1 * l2 * x(6) * sin(x(3) - x(2)), 0];
 
g = [0; m12 * g * l1 * sin(x(2)); m2 * g * l2 * sin(x(3))];

h = [1 0 0]';

syms u

f = [zeros(3), eye(3); zeros(3), -inv(D)*C] * x ...
    + [zeros(3, 1); -inv(D) * g] + [zeros(3, 1);-inv(D)*h]* u;

%% Animation
