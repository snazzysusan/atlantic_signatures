function dx=atlantic_function_mod(T,Y)
%This function is called by the my_phase_plane function to be solved using 
%ode45 and plotted as x vs. xdot on a phase plane.
x=Y(1); % Variable vectors
y=Y(2);

%ocean current
    fluid_rot = 5*pi/180;
    x_fac = 2;
    y_fac = 1;
    v_theta = -2;

    %to clean up the code
    k1 = x_fac.*cos(fluid_rot);
    k2 = y_fac.*sin(fluid_rot);
    k3 = x_fac.*sin(fluid_rot);
    k4 = y_fac.*cos(fluid_rot);

    den = sqrt(x^2+y^2);

%Magnetic Field
    %parameters
    gamma_rot_angle = 10;
    beta_gamma_space = 5;
    beta_rot_angle = -(90 - gamma_rot_angle - beta_gamma_space);
    gamma_rot_angle = gamma_rot_angle*pi/180;
    beta_rot_angle = beta_rot_angle*pi/180;

    beta0 = [0;0;0];
    gamma0 = [0;0;0];

    beta_normal = [1; 0; 1];
    gamma_normal = [0; -1; 1];

    regular_sensing = 1;
    x_initial = [-20, 5];

    distance_beta = sum(beta_normal.*beta0);
    distance_gamma = sum(gamma_normal.*gamma0);
    
    %Magnetic Field Functions
    gamma = @(X,Y, lambda) (distance_gamma - gamma_normal(1).*(X*cos(lambda) + Y.*sin(lambda)) -...
        gamma_normal(2).*(-X*sin(lambda) + Y.*cos(lambda)))./gamma_normal(3);
    beta = @(X,Y,lambda) regular_sensing*(distance_beta - beta_normal(1).*(X*cos(lambda) + Y.*sin(lambda)) -...
        beta_normal(2).*(-X*sin(lambda) + Y.*cos(lambda)))./beta_normal(3);

    gamma_plot = gamma(x,y,gamma_rot_angle);
    beta_plot = beta(x,y,beta_rot_angle);

    beta_goal = beta(x_initial(:,1), x_initial(:,2), beta_rot_angle);
    gamma_goal = gamma(x_initial(:,1), x_initial(:,2), gamma_rot_angle);
    
    %used in the final system (dx)
    d_beta = beta_goal - beta_plot;
    d_gamma = gamma_goal - gamma_plot;
    d_magnitude = sqrt((d_beta).^2 + (d_gamma).^2);


%Full System
%Using cartesian instead of polar speeds up the ODE Solver a lot
dx=zeros(2,1);
dx(1)=((d_beta/d_magnitude)) - (v_theta/den)*(k1*y + k2*x);
dx(2)=((d_gamma/d_magnitude)) + (v_theta/den)*(k4*x - k3*y);

%Original Code
%dx=zeros(2,1)
% dx(1)=((1/d_magnitude)*d_beta)-x_fac*sin(atan(y/x))*v_theta*cos(fluid_rot)-(y_fac*cos(atan(y/x))...
%     *sin(fluid_rot)*v_theta);
% dx(2)=((1/d_magnitude)*d_gamma)-x_fac*sin(atan(y/x))*v_theta*sin(fluid_rot)+(y_fac*cos(atan(y/x))...
%     *cos(fluid_rot)*v_theta);


end