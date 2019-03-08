function [ zk, Hk ] = insErrorMeasurementModel_GNSS_rv( xk, wk )
% Measurement model for GNSS Position and Velocity error components
% Computes the measurement model of the errors for the INS - GPS state 
% difference.
% 
% This models the 21 element error state of the INS, includes modeling
% the accel and gyro scale factor errors.
% 
% INPUTS
% xk - 21 x 1 double vector
%      State vector at time k
%          [ del_r_e;       - Position Error
%            del_v_e;       - Velocity Error
%            err_e;         - Misalignment Error
%            del_biasacc_e; - Acc Bias Drift Error
%            del_biasgyr_e; - Byr Bias Drift Error
%            biasinitacc_e; - Initial Acc Bias Error
%            biasinitgyr_e; - Initial Gyr Bias Error
%            scaleacc_e;    - Acc Scale Factor Error
%            scalegyr_e ]   - Gyr Scale Factor Error
% wk - 6 x 1 double vector
%      Process Noise vector at time k
%          [ nu_r_GNSS_x;  - GNSS position measurement noise
%            nu_r_GNSS_y;  - GNSS position measurement noise
%            nu_r_GNSS_z;  - GNSS position measurement noise
%            nu_v_GNSS_x;  - GNSS velocity measurement noise
%            nu_v_GNSS_y;  - GNSS velocity measurement noise
%            nu_v_GNSS_z ] - GNSS velocity measurement noise
% 
% OUTPUTS
% zk - 6 x 1 double vector
%      State vector at time k
% Hk - 6 x 21 double matrix
%      State Transition partial derivative
% 
% @author: Matt Marti
% @date: 2019-03-08

nx = 21;
nz = 6;

Hk = zero(nz,nx);
Hk(1:nz,1:nz) = eye(nz);

zk = Hk * xk + wk;

end

