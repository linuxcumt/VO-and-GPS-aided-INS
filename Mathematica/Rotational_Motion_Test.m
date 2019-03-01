
function temp

    % Parameters
    dx10 = .2;
    ddxpsi = 0.4;
    psi0 = 1;
    psidot = 2;

    t = 0:0.1:2;
    
    % Sine (x1)
    f1 = @(t) dx10 + 1/psidot * ddxpsi .* sin(psi0 + psidot * t);
%     f2 = @(t) dx10 + ddxpsi.*t.*sa(0.5*psidot*t).*cos(psi0 + 0.5*psidot*t);
    f3 = @(t) dx10 + ddxpsi.*(sin(psi0)*cos(psidot*t)/psidot + cos(psi0)*sin(psidot*t)/psidot);
        
    % Velocity why is it wrong
    figure(1)
    hold off
%     plot(t,[f1(t);f2(t);f3(t)], 'linewidth', 2)
    plot(t,[f1(t);f3(t)], 'linewidth', 2)
    title('x_1 Velocity')
    
    figure(2)
    hold off
%     plot(t,[f1(t)-f2(t);f1(t)-f3(t)], 'linewidth', 2)
    plot(t,f1(t)-f3(t), 'linewidth', 2)
    title('x_1 Velocity Error')
    
    % Velocity True
    f2 = @(t) dx10 + ddxpsi.*t.*sa(0.5*psidot*t).*cos(psi0 + 0.5*psidot*t);
    f4func = @(t) ddxpsi.*cos(psi0+psidot*t);
    f4 = @(t) int(f4func,t,0) + dx10;
    
    figure(3)
    hold off
    plot(t,[f2(t);f4(t)], 'linewidth', 2)
    title('x_1 Velocity')
    
    % Position
    f5 = @(t) dx10*t - ddxpsi.*(cos(psi0+0.5*psidot*t)+t.*psidot.*sin(psi0))/(psidot^2);
    f6func = @(t) dx10 + ddxpsi.*t.*sa(0.5*psidot*t).*cos(psi0 + 0.5*psidot*t);
    f6 = @(t) int(f6func,t,0);
    
    figure(4)
    hold off
    plot(t,[f5(t);f6(t)], 'linewidth', 2)
    title('x_1 Position')
end

function y = int(f,t,y0)

y = zeros(size(t));
y(1) = y0;
for i = 2:length(t)
    t1 = t(i-1);
    t2 = t(i);
    f1 = f(t1);
    f2 = f(t2);
    
    dy = f1*(t2-t1) + 0.5*(f2-f1)*(t2-t1);
    y(i) = y(i-1) + dy;
end

end

function r = sa(x)
r = x;
for i = 1:length(x)
    if ~x(i)
        r(i) = 1;
    else
        r(i) = sin(x(i))/x(i);
    end
end
end

function r = ca(x)
r = x;
for i = 1:length(x)
    if ~x(i)
        r(i) = 1;
    else
        r(i) = cos(x(i))/x(i);
    end
end
end