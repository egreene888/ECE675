%% Funwork 2
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

%% Linearizing the model around a non-zero equilibrium 
%
% We can use the Taylor Linearization to approximate the equations of 
% motion about a point $x_e, \, u_e$. 
% 
% $$ 
% \mathbf{\dot x} \approx \mathbf f(\mathbf x_e, u_e)
% + \frac{\partial \mathbf f}{\partial \mathbf x} (\mathbf x - \mathbf x_e)
% + \frac{\partial \mathbf f}{\partial u} (u - u_e)
% $$
