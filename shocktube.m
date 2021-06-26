% Note: No EOS for saturated magma is specified. As such, the saturated
% magma in Turcotte et al. is really a magma-vapour mixture, in the limit
% of pressure p -> p_0 (the saturation pressure) such that all the
% volatiles are dissolved in the melt.
% For Wilson&Head values for H2O in basaltic magma, the saturation pressure
% is typically below atmospheric. If we want a physically sound problem, we
% need to supply low enough pressures that (only) the vapour phase exists,
% or specify a large mass concentration. For p0 = 1 atm, we need X0 = 15.2
% which means volatile mass is 15 times that of the melt.
%
% No EOS is specified for the magma, so the fluid pressure is unknown as a
% function of density, (and possibly temperature). In the Turcotte model
% this is worked around by having only saturated magma at pressure exactly
% p0, so the magma is simply the upper-bound pressure of the magma-vapour
% mixture.

%% Physics parameters
% Isothermal parameters
R = 8.314/44.01*1e3;    % CO2
T0 = 1500;  % K
E0 = R * T0;
% Magma mixture properties
rho10 = 2.7; % Dry magma density (kept low for a smaller contrast)
X0 = 12;  % Mass concentration of volatile
% Henry's law for H2O in basaltic magma Wilson & Head (1981)
k = 2.3e-6;
n = 1;
p0 = (X0/k)^(1.0/n);
%% Meshing
Ncores = 12;
mesh.x = linspace(-100,100,12*Ncores)';
dx = mesh.x(2) - mesh.x(1);
dt = dx/8000;
tFinal = 5e-2;
tVec = 0:dt:tFinal;
%% Set initial conditions
f0 = (1-1e-5)*(mesh.x < 0);
u0 = zeros(size(mesh.x));
pg = 101.3e3;          % Yes, the gas pressure (pg) must be quite small
rhog0 = pg / (R * T0);   
rho0 = rhog0 * ones(size(mesh.x)); % Gas density fill
rho0(mesh.x < 0) = rho10*(1 + X0);     % Replace with mixture density
q0 = [rho0, rho0.*u0, f0];
%% Forward Euler
methods = {'MatrixSplit-0', 'LF-1', 'LLF-0'};
soln = {};
timings = {};
figure(777); clf;
for idxMethod = 1:length(methods)
    timerStart = tic;
    method = methods{idxMethod};
    % Initialize data storage
    qData = nan(length(mesh.x),3,length(tVec));
    qData(:,:,1) = q0;
    % Initialize solution
    q = q0;
    for i = 2:length(tVec)
        t = tVec(i);
        % Forward Euler
        [fluxes, maxEig] = computeFluxes(dx, q, E0, rho10, X0, n, p0, ...
            dt, method);
        if maxEig > dx/dt
            error("CFL condition violated: c dt / dx == " ...
                  + maxEig * dt / dx)
        end
        % Step
        q = q - dt/dx*fluxes;
        % Save data
        qData(:,:,i) = q;
        % Plot states
        clf;
        plotState(mesh.x, q, E0, rho10, X0, n, p0);
        drawnow;
        title("t = " + t);
    end
    soln{idxMethod} = qData;
    timings{idxMethod} = toc(timerStart);
end

%% Post-process plotting
figure(778); clf;
plot(mesh.x, soln{1}(:,1,end), '.-', 'LineWidth', 1); hold on
plot(mesh.x, soln{2}(:,1,end), '.-', 'LineWidth', 1);
plot(mesh.x, soln{3}(:,1,end), '.-', 'LineWidth', 1); hold off
xlabel 'x [m]'
ylabel '\rho [kg/m^3]'
title("t = " + tFinal)
legend({'Matrix Split p0', 'Lax-Friedrichs p1 (Van Leer)', 'Rusanov p0'})
grid on
grid minor

%% View density
rhoFinal = soln{1}(:,1,end);
uFinal = soln{1}(:,2,end)./soln{1}(:,1,end);
fFinal = soln{1}(:,3,end);
pFinal = arrayfun(@(rho, u, f) pFn(rho, u, f, E0, rho10, X0, n, p0), ...
        q(:,1), q(:,2)./q(:,1), q(:,3) );
    
figure(779); clf;
subplot(4,1,1);
plot(mesh.x, rhoFinal, '.-', 'LineWidth', 1);
xlabel 'x [m]'
ylabel '\rho [kg/m^3]'
title("t = " + tFinal + " s")
grid on
grid minor

subplot(4,1,2);
plot(mesh.x, uFinal, '.-', 'LineWidth', 1);
xlabel 'x [m]'
ylabel 'u [m/s]'
title("t = " + tFinal + " s")
grid on
grid minor

subplot(4,1,3);
plot(mesh.x, fFinal, '.-', 'LineWidth', 1);
xlabel 'x [m]'
ylabel 'f'
title("t = " + tFinal + " s")
grid on
grid minor

subplot(4,1,4);
plot(mesh.x, pFinal, '.-', 'LineWidth', 1);
xlabel 'x [m]'
ylabel 'p'
title("t = " + tFinal + " s")
grid on
grid minor

%% Nondimensional plots
nondim.u = uFinal * sqrt(rho10 / p0);
nondim.eps = p0 / rho10 / R / T0;
nondim.p = pFinal / p0;

% Get rarefaction solution from theory
theory.x = sqrt(p0/rho10) * tFinal * ...
    (-coeff*log(nondim.p) - ...
    (X0 + nondim.p*(nondim.eps - X0))/sqrt(nondim.eps*X0*(1+X0)));
% Fill with unperturbed magma to the left of rarefaction fan extent
theory.x = [mesh.x(1); theory.x];
% Reuse p from numeric solution only as the range
theory.p = [1; nondim.p];

figure(780); clf;
subplot(1,2,1);
coeff = sqrt(X0 / nondim.eps / (1 + X0));
plot(mesh.x, nondim.u, '.-', 'LineWidth', 1);
hold on
plot(mesh.x, -coeff*log(nondim.p), '.-', 'LineWidth', 1)
legend({'$u / \sqrt{\rho_{10} / p_0}$', ...
    '$-[ \frac{X_0}{\varepsilon(1 + X_0)} ]^{1/2} \ln \frac{p}{p_0}$'}, ...
    'Interpreter', 'latex', 'FontSize', 14, 'location', 'northwest')
xlabel ('$x$ [m]', 'Interpreter', 'latex')
ylabel ('$u / \sqrt{\rho_{10} / p_0}$', 'Interpreter', 'latex')
title ("t = " + tFinal +" s")

subplot(1,2,2);
plot(mesh.x, nondim.p, '.-', 'LineWidth', 1);
hold on
plot(theory.x, theory.p, '-', 'LineWidth', 1);
legend({'Numerical solution', ...
    'Exact expansion fan, $n = 1$'}, ...
    'Interpreter', 'latex', 'FontSize', 12, 'location', 'northeast')
xlabel ('$x$ [m]', 'Interpreter', 'latex')
ylabel ('$p / p_0$', 'Interpreter', 'latex')
title ("t = " + tFinal +" s")

%% Examine flux Jacobian eigenvalues
% Extra post process step to look at eigenvalues at some point in the final
% solution.
qQuery = soln{1}(length(mesh.x)/2,:,end);
rhoQuery = qQuery(1);
uQuery = qQuery(2)/rhoQuery;
fQuery = qQuery(3);

A = fluxJacobian(qQuery(1), qQuery(2)/qQuery(1), qQuery(3), ...
    E0, rho10, X0, n, p0);
[eigVecs, eigVals] = eig(A);

%% Plotting function
function plotState(x, q, E0, rho10, X0, n, p0)
    subplot(2,2,1);
    plot(x, q(:,1), '.-');
    xlabel 'x'
    ylabel '\rho'
    yL = ylim;
    ylim([0, yL(2)])

    subplot(2,2,2);
    plot(x, q(:,2) ./ q(:,1), '.-');
    xlabel 'x'
    ylabel 'u'

    subplot(2,2,3);
    plot(x, q(:,3), '.-');
    xlabel 'x'
    ylabel 'f'
    
    subplot(2,2,4);
    plot(x, arrayfun(@(rho, u, f) pFn(rho, u, f, E0, rho10, X0, n, p0), ...
        q(:,1), q(:,2)./q(:,1), q(:,3) ), '.-');
    xlabel 'x'
    ylabel 'p'
end

%% Pressure computation [implicit, nonlinear solve]
function p_ = pFn(rho, u, f, E0, rho10, X0, n, p0)
    if f >= 1
        p_ = p0;
        return
    end
    % Computes pressure
    pGuess = p0;
    pFn = @(pVar) pVar - E0 ./ (1.0 - f) .* ...
        (rho - rho10 * f * (1 + X0 * (pVar./p0).^n));
%     searchAxis = linspace(0,1000,1000);
%     plot(searchAxis, pFn(searchAxis));
    [p_, ~, exitflag] = fsolve( ...
        pFn,pGuess,optimoptions('fsolve','Display','off', ...
            'MaxIterations', 100000));
    if exitflag ~= 1 && exitflag ~= 2
        % Non-physical state (probably p < 0, meaning more exsolution
        % than possible). Consider using lower contrast in the initial.
        error("Extrapolated to a negative pressure." + ...
              "More exsolution than possible." + ...
              "Initial condition may be unsound.")
    end
    if p_ > p0
        % This may be required briefly for the initial condition f = 1-eps
        % because the pressure is highly sensitive to f near 1
        % error('Unhandled pressure case.')
    end
end

function A = fluxJacobian(rho, u, f, E0, rho10, X0, n, p0)
% Computes the 3x3 flux Jacobian.
function y = denom(rho, u, f, E0, rho10, X0, n, p0)
    y = 1.0 + E0 * f ./ (1 - f) .* rho10 .* X0 .* n .* ...
        (pFn(rho, u, f, E0, rho10, X0, n, p0)./p0).^n ./ ...
        pFn(rho, u, f, E0, rho10, X0, n, p0);
end
function y = dpdrho(rho, u, f, E0, rho10, X0, n, p0)
    y = 1 ./ denom(rho, u, f, E0, rho10, X0, n, p0) .* ...
        E0 ./ (1.0 - f);
end
function y = dpdf(rho, u, f, E0, rho10, X0, n, p0)
    y = 1 ./ denom(rho, u, f, E0, rho10, X0, n, p0) .* ...
        (E0 ./ (1 - f).^2) .* (rho - rho10 * (1 + X0 * ...
        (pFn(rho, u, f, E0, rho10, X0, n, p0)./p0).^n));
end
A = [0, 1, 0;
     -u.^2 + dpdrho(rho, u, f, E0, rho10, X0, n, p0), 2*u, ...
         dpdf(rho, u, f, E0, rho10, X0, n, p0);
     -f .* u ./ rho, f ./ rho, u];

% % Wavespeed
% c = sqrt(f / rho * dpdf(rho, u, f, E0, rho10, X0, n, p0) + ...
%      dpdrho(rho, u, f, E0, rho10, X0, n, p0));
end

function fl = fluxFn(rho, u, f, E0, rho10, X0, n, p0)
% Compute 3x1 flux vector
    fl = [rho .* u;
          rho .* u.^2 + pFn(rho, u, f, E0, rho10, X0, n, p0);
          f .* u];
end

function [dflux, eigLargest] = computeFluxes(dx, q, E0, rho10, X0, n, ...
                                             p0, dt, method)
    % Compute numerical flux F (3x1 vector).
    % Assumes Dirichlet boundaries (don't let the waves hit the boundary!)

    % Initialize
    dflux = zeros(size(q));
    eigLargest = 0;
    
    % Skip first-order reconstruction for select methods
    if strcmpi('MatrixSplit-0', method) ...
            || strcmpi('LLF-0', method)
        skipRecon = true;
    else
        skipRecon = false;
    end
    
    %% Slope limiters
    function phi = limiterVL(r)
        if isinf(r)
            phi = 2;
        elseif isnan(r)
            phi = 2; % Constant case
        elseif r <= 0
            phi = 0;
        else
            phi = (r+abs(r))/(1+r);
        end
    end

    function phi = limiterPlus(uLeft, uCenter, uRight)
        r = (uRight - uCenter)/(uCenter - uLeft);
        phi = limiterVL(r);
    end

    function phi = limiterMinus(uLeft, uCenter, uRight)
        r = (uCenter - uLeft)/(uRight - uCenter);
        phi = limiterVL(r);
    end

    if ~skipRecon
        %% Compute slope-limited piecewise linear reconstruction at edges
        % Compute state vectors of size N-2, skipping the boundary cells
        % and the first cell on the left for leftward bias (resp. the first
        % cell on the right for rightward bias).
        
        rho = q(:,1);
        u = q(:,2) ./ q(:,1);
        f = q(:,3);
        
        for i = 2:size(q,1)-1
            %% Compute left-biased edge states W_{i+1/2}
            % Limiter by Van Leer
            rhoEdge_LeftBiased(i) = rho(i) + 0.5 * ...
                limiterPlus(rho(i-1), rho(i), rho(i+1))*(rho(i)-rho(i-1));
            uEdge_LeftBiased(i) = u(i) + 0.5 * ...
                limiterPlus(u(i-1), u(i), u(i+1))*(u(i)-u(i-1));
            fEdge_LeftBiased(i) = f(i) + 0.5 * ...
                limiterPlus(f(i-1), f(i), f(i+1))*(f(i)-f(i-1));
            %% Compute right-biased edge states W_{i-1/2}
            rhoEdge_RightBiased(i) = rho(i) - 0.5 * ...
                limiterMinus(rho(i-1), rho(i), rho(i+1))*(rho(i+1)-rho(i));
            uEdge_RightBiased(i) = u(i) - 0.5 * ...
                limiterMinus(u(i-1), u(i), u(i+1))*(u(i+1)-u(i));
            fEdge_RightBiased(i) = f(i) - 0.5 * ...
                limiterMinus(f(i-1), f(i), f(i+1))*(f(i+1)-f(i));
        end
        %% Extend edge states to boundaries assuming constant-flux BCs
        % Take N-2 biased edge states to N+1 edge states
        extendDataL = @(vec) [vec(2), vec(2), vec(2:end), vec(end)]';
        extendDataR = @(vec) [vec(2), vec(2:end), vec(end), vec(end)]';
        rhoEdge_LeftBiased = extendDataL(rhoEdge_LeftBiased);
        uEdge_LeftBiased = extendDataL(uEdge_LeftBiased);
        fEdge_LeftBiased = extendDataL(fEdge_LeftBiased);
        rhoEdge_RightBiased = extendDataR(rhoEdge_RightBiased);
        uEdge_RightBiased = extendDataR(uEdge_RightBiased);
        fEdge_RightBiased = extendDataR(fEdge_RightBiased);
        % Compute conservative variables
        Q_LeftBiased = [rhoEdge_LeftBiased, ...
                        rhoEdge_LeftBiased .* uEdge_LeftBiased, ...
                        fEdge_LeftBiased];
        Q_RightBiased = [rhoEdge_RightBiased, ...
                         rhoEdge_RightBiased .* uEdge_RightBiased, ...
                         fEdge_RightBiased];
    end
    
    if strcmpi('MatrixSplit-0', method)
        % Steger-Warming-type matrix split upwinding with piecewise
        % constant reconstruction
        parfor i = 2:size(q,1)-1
            qLoc = q(i,:)';
            rho = qLoc(1);
            u = qLoc(2) / qLoc(1);
            f = qLoc(3);

            frozenCoeff = false;
            if frozenCoeff
                % Frozen coefficient problem:
                A = fluxJacobian(1.2, 0, 0, E0, rho10, X0, n, p0);
                [R, L] = eig(A);
                L(3,3) = 0.5*(abs(L(1,1)) + abs(L(2,2)));
                A = R * L / R;
    %             A = eye(3,3);
            else
                A = fluxJacobian(rho, u, f, E0, rho10, X0, n, p0);
            end

            % Naive A matrix-split
            [R, L] = eig(A);
            eigLargest = max(eigLargest, max(abs(diag(L))));
            % Harten Entropy-fix parameter
            reg = 1e-8;
            % Extract positive and negative parts of the matrix
            Lplus = 0.5 * (L + sqrt(L.^2 + reg.^2*eye(3)));
            Lminus = 0.5 * (L - sqrt(L.^2 + reg.^2*eye(3)));
    %         Aplus = 0.5 * R * (L + sqrt(L.^2 + reg.^2*eye(3))) / R;
    %         Aminus = 0.5 * R * (L - sqrt(L.^2 + reg.^2*eye(3))) / R;

            %% Gradient computation
            gradLeft = (qLoc - q(i-1,:)')/dx;
            gradRight = (q(i+1,:)' - qLoc)/dx;
            FPlusLoc = R * (Lplus * (R \ gradLeft));
            FMinusLoc = R * (Lminus * (R \ gradRight));
            dflux(i,:) = dx*(FPlusLoc + FMinusLoc);
        end
    elseif strcmpi('LF-1', method) % 
        % Lax-Friedrichs two-wave approx with piecewise linear
        % reconstruction
        parfor i = 2:size(q,1)-1
            %% Flux function
            fluxCompute = @(q) fluxFn(q(1), q(2)/q(1), q(3), E0, rho10, ...
                X0, n, p0);
            % Left flux
            qLeftBiased = Q_LeftBiased(i,:)';
            qRightBiased = Q_RightBiased(i,:)';
            fL = fluxCompute(qLeftBiased);
            fR = fluxCompute(qRightBiased);
            % Lax Friedrichs flux
            fluxLeft = 0.5*(fL + fR) - 0.5*dx/dt *...
                (qRightBiased - qLeftBiased);
            % Right flux
            qLeftBiased = Q_LeftBiased(i+1,:)';
            qRightBiased = Q_RightBiased(i+1,:)';
            fL = fluxCompute(qLeftBiased);
            fR = fluxCompute(qRightBiased);
            % Lax Friedrichs flux
            fluxRight = 0.5*(fL + fR) - 0.5*dx/dt *...
                (qRightBiased - qLeftBiased);
            % Compute total flux
            dflux(i,:) = fluxRight - fluxLeft;
        end
    elseif strcmpi('LLF-0', method) 
        % Local Lax-Friedrichs (Rusanov) with piecewise constant
        % reconstruction
        parfor i = 2:size(q,1)-1
            %% Largest eigenvalue function
            eigLargestFn = @(q) max(abs(eig(fluxJacobian(...
                q(1), q(2)/q(1), q(3), ...
                E0, rho10, X0, n, p0))));
            %% Flux function
            fluxCompute = @(q) fluxFn(q(1), q(2)/q(1), q(3), E0, rho10, ...
                X0, n, p0);
            % Left flux
%             qLeftBiased = Q_LeftBiased(i,:)';
%             qRightBiased = Q_RightBiased(i,:)';
%             fL = fluxCompute(qLeftBiased);
%             fR = fluxCompute(qRightBiased);
            fL = fluxCompute(q(i-1,:));
            fR = fluxCompute(q(i,:));
            % Local speed parameter (Rusanov)
%             c = max([eigLargestFn(q(i-1,:)), ...
%                      eigLargestFn(q(i,:)), ...
%                      eigLargestFn(Q_LeftBiased(i,:)), ...
%                      eigLargestFn(Q_RightBiased(i,:))]);
            c = max([eigLargestFn(q(i-1,:)), ...
                     eigLargestFn(q(i,:))]);
            % Lax Friedrichs flux
%             fluxLeft = 0.5*(fL + fR) - 0.5*c *(qRightBiased - qLeftBiased);
            fluxLeft = 0.5*(fL + fR) - 0.5*c *(q(i,:) - q(i-1,:))';
            % Right flux
%             qLeftBiased = Q_LeftBiased(i+1,:)';
%             qRightBiased = Q_RightBiased(i+1,:)';
%             fL = fluxCompute(qLeftBiased);
%             fR = fluxCompute(qRightBiased);
            fL = fR; %fluxCompute(q(i,:));
            fR = fluxCompute(q(i+1,:));
            % Local speed parameter
%             c = max([eigLargestFn(q(i,:)), ...
%                      eigLargestFn(q(i+1,:)), ...
%                      eigLargestFn(Q_LeftBiased(i+1,:)), ...
%                      eigLargestFn(Q_RightBiased(i+1,:))]);
            c = max([eigLargestFn(q(i,:)), ...
                     eigLargestFn(q(i+1,:))]);
            % Lax Friedrichs flux
            fluxRight = 0.5*(fL + fR) - 0.5*c *(q(i+1,:) - q(i,:))';
            % Compute total flux
            dflux(i,:) = fluxRight - fluxLeft;
        end
    else
        error('No flux')
    end
    
end