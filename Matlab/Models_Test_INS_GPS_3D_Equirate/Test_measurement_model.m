%% Test_measurement_model.m
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
nz = 6;

for ntests = 1:100
    
    % State
    xk = 1000 * (rand(nx,1) - 0.5);
    wk = 1000 * (rand(nz,1) - 0.5);

    % Numerical derivative H
    dx = 1e-4;
    H_num = zeros(nz,nx);
    for i = 1:nx
        xp = xk;
        xm = xk;
        xp(i) = xp(i) + dx;
        xm(i) = xm(i) - dx;
        zpkp1 = measurement_model( xp, wk );
        zmkp1 = measurement_model( xm, wk );
        H_num(:,i) = (zpkp1 - zmkp1) / (2*dx);
    end

    % Measurement Model
    [z, H] = measurement_model( xk, wk );

    % Assertions
    p = 2e-6;
    H_diff = H - H_num;
    Hcheck = abs(H_diff) <= p;
    assert(sum(sum(Hcheck)) == nx*nz, 'Bad H');
end


%% Output
fprintf('PASSED: Test_measurement_model\n');