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
% where $K_0$ is the kinetic energy due to the lateral motion of the cart 
% and pendulum, $K_1$ is due to the rotation of mass 1 about its pivot point
% $p$, and $K_2$ is due to the rotation of mass 2 about $p$. The value of 
% $K_0$ is 
% 
% $$ K_0 = (M+m_1+m_2)\ddot{x} $$
% 
% While the value of $K_1$ is 
% 
% $$
% \begin{array}{rcl}
% K_1 & = & \frac{1}{2}m_1 v_1^2 \\
% & = & m_1 l_1^2 \dot{\theta_1}^2
% \end{array}
% $$
% 
% and the value of $k_2$ is 
% 
% $$ 
% \begin{array}{rcl}
% K_2 & = & \frac{1}{2}m_2 v_2^2 \\
% & = & \frac{1}{2} m_2 \left[ (-l_1 \dot{\theta_1} \sin \theta_1 - 
% l_2 \dot{\theta}_2 \sin \theta_2)^2 
% + (l_1 \dot{\theta_1} \cos \theta_1 + 
% l_2 \dot{\theta}_2 \cos \theta_2)^2 \right] \\
% & = & \frac{1}{2}m_2 \left[ l_1^2 \dot{\theta}_1^2 + 
% l_2^2 \dot{\theta}_2^2 + 2 l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 
% \cos(\theta_2 - \theta_1) \right] \\
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
% L = \frac{1}{2} \left(M + m_1 + m_2\right) \dot{x}^2 + 
% \frac{1}{2} m_1 l_1^2 \dot{\theta}_1^2  +
% \frac{1}{2} m_2 \left[ l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +
% 2 l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos(\theta_2 - \theta_1) \right]
% - m_1 g l_1 \cos \theta_1 - 
% m_2 g \left( l_1 \cos \theta_1 + l_2 \cos \theta_2 \right)
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
% \frac{\partial L}{\partial \dot{x}} = (M + m_1 + m_2)\dot{x} 
% $$
% 
% $$
% \frac{d}{dt}\frac{\partial L}{\partial \dot{x}} = (M + m_1 + m_2)\ddot{x}
% $$
% 
% $$ 
% \frac{\partial L}{\partial x} = 0
% $$
% 
% Therefore, the first equation of motion is 
% 
% $$ (M + m_1 + m_2)\ddot{x} = u $$
% 
% We can find the second equation of motion by finding
% 
% $$
% \frac{\partial L}{\partial \dot{\theta}_1} = 
% (m_1 + m_2) l_1^2\dot{\theta}_1 + 
% m_2 l_1 l_2 \dot{\theta}_2 \cos(\theta_2 - \theta_1)
% $$
% 
% $$ 
% \begin{array}{rcl} 
% \frac{d}{dt} \frac{\partial L}{\partial \dot{\theta}_1} & = & 
% (m_1 + m_2) l_1^2 \ddot{\theta_1} + 
% m_2 l_1 l_2 \ddot{\theta}_2 \cos (\theta_2 - \theta_1) - 
% m_2 l_1 l_2 \dot{\theta}_2 (\dot{\theta}_2 - 
% \dot{\theta}_1) \sin (\theta_2 - \theta_1) \\
% & = & (m_1 + m_2) l_1^2 \ddot{\theta}_1 + 
% m_2 l_1 l_2 \ddot{\theta}_2 \cos (\theta_2 - \theta_1) - 
% m_2 l_1 l_2 \dot{\theta}_2^2 \sin (\theta_2 - \theta_1) + 
% m_2 l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \sin (\theta_2 - \theta_1) 
% \end{array}
% $$