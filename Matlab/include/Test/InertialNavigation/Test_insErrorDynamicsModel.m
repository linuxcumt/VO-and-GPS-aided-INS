%% Test_insErrorDynamicsModel
%
% Tests the 27 state error propagation model
%
% @author: Matt Marti
% @date: 2019-03-09

clear; clear global
constants;


%% Test

for n = 1:100

    % Parameters
    dt = 0.01;
    beta_acc = 0.01*rand(3,3);
    beta_gyr = 0.01*rand(3,3);
    beta_sf_acc = 0.01*rand(3,3);
    beta_sf_gyr = 0.01*rand(3,3);
    nx = 27;
    nv = 18;
    nInt = 100;

    % INS State
    ins = zeros(22,1);
    ins(1:3,1) = 6700e3*[1; 1; 0]; % Arbitrary position vector
    ins(4:6,1) = [2; 0; 0]; % Arbitrary velocity vector
    q = [1; 0; 0; 0]; % Set heading, pitch, roll to 0
    q = q/norm(q); % Normalize quaternion
    ins(7:10,1) = q; % Attitude

    % Error Estimates (the state)
    xk = zeros(nx,1);
    xk(1:3,1) = 5*randn(3,1); % Position error
    xk(4:6,1) = 0.1*randn(3,1); % Velocity error
    xk(7:9,1) = 0.01*randn(3,1); % Tilt error
    xk(10:12,1) = 0.01*randn(3,1); % Bias Drift Rate acc error
    xk(13:15,1) = 0.01*randn(3,1); % Bias Drift Rate gyr error
    xk(16:18,1) = 0.01*randn(3,1); % Bias acc error
    xk(19:21,1) = 0.01*randn(3,1); % Bias gyr error
    xk(22:24,1) = 0.01*randn(3,1); % Scale Factor acc error
    xk(25:27,1) = 0.01*randn(3,1); % Scale Factor gyr error
    
    % Noise vector
    vk = zeros(nv,1);

    % Run function
    [ xkp1, Fk, Gammak ] ...
        = insErrorDynamicsModel( dt, xk, vk, ins, beta_acc, beta_gyr, ...
                                 beta_sf_acc, beta_sf_gyr, nInt );

    % Numerical derivative - Fk
    dx = 1e-4;
    D_num = zeros(nx,nx);
    for i = 1:nx
        xp = xk;
        xm = xk;
        xp(i) = xp(i) + dx;
        xm(i) = xm(i) - dx;
        xkp1p = insErrorDynamicsModel( dt, xp, vk, ins, beta_acc, beta_gyr, ...
                                       beta_sf_acc, beta_sf_gyr, nInt );
        xkp1m = insErrorDynamicsModel( dt, xm, vk, ins, beta_acc, beta_gyr, ...
                                       beta_sf_acc, beta_sf_gyr, nInt );
        D_num(:,i) = (xkp1p - xkp1m) ./ (2*dx);
    end

    % Assertions - Fk
    p = 2e-6;
    D_diff = Fk - D_num;
    Dcheck = abs(D_diff) <= p;
    assert(prod(prod(Dcheck)) == 1, 'Bad Fk Derivative');

    % Numerical derivative - Gammak
    dv = 1e-4;
    D_num = zeros(nx,nv);
    for i = 1:nv
        vp = vk;
        vm = vk;
        vp(i) = vp(i) + dv;
        vm(i) = vm(i) - dv;
        xkp1p = insErrorDynamicsModel( dt, xk, vp, ins, beta_acc, beta_gyr, ...
                                       beta_sf_acc, beta_sf_gyr, nInt );
        xkp1m = insErrorDynamicsModel( dt, xk, vm, ins, beta_acc, beta_gyr, ...
                                       beta_sf_acc, beta_sf_gyr, nInt );
        D_num(:,i) = (xkp1p - xkp1m) ./ (2*dx);
    end

    % Assertions - Gammak
    p = 2e-6;
    D_diff = Gammak - D_num;
    Dcheck = abs(D_diff) <= p;
    assert(prod(prod(Dcheck)) == 1, 'Bad Gammak Derivative');
end


%% Output
fprintf('PASSED: Test_insErrorDynamicsModel\n');