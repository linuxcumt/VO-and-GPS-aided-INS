%% Test_dynamics_model.m
% 
% Script to test partial derivatives of the INS Dynamics Model
% function. This was developed for the VO and GPS aided INS project in the
% VT AOE 5984 Machine Learining class.
% 
% @author: Matt Marti
% @date: 2019-03-01

clear


%% Test

% Dynamics model
nx = 21;
nv = 12;

% Do numerical tests 100 times
for ntests = 1:100
    
    % State
    xk = 1000 * (rand(nx,1) - 0.5);
    vk = 100 * (rand(nv,1) - 0.5);

    % Run function for exact values
    [ xkp1, F, Gamma ] = dynamics_model( xk, vk );

    % Numerical derivative F
    dx = 1e-4;
    F_num = zeros(nx,nx);
    for i = 1:nx
        xp = xk;
        xm = xk;
        xp(i) = xp(i) + dx;
        xm(i) = xm(i) - dx;
        xpkp1 = dynamics_model( xp, vk );
        xmkp1 = dynamics_model( xm, vk );
        F_num(:,i) = (xpkp1 - xmkp1) / (2*dx);
    end

    % Numerical derivative Gamma
    dv = dx;
    Gamma_num = zeros(nx,nv);
    for i = 1:nv
        vp = vk;
        vm = vk;
        vp(i) = vp(i) + dx;
        vm(i) = vm(i) - dx;
        xpkp1 = dynamics_model( xk, vp );
        xmkp1 = dynamics_model( xk, vm );
        Gamma_num(:,i) = (xpkp1 - xmkp1) / (2*dx);
    end

    % Assertions
    p = 1e-9;
    F_diff = F - F_num;
    Gamma_diff = Gamma - Gamma_num;
    Fcheck = abs(F_diff) <= p;
    assert(sum(sum(Fcheck)) == nx*nx, 'Bad F');
    Gammacheck = abs(Gamma_diff) <= p;
    assert(sum(sum(Gammacheck)) == nv*nx, 'Bad Gamma');
end

%% Output
fprintf('PASSED: Test_dynamics_model\n');