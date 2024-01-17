function Base_Flow = get_base_flow(Definitions)

% Get data from user input / fixed user input
lambdaGuesses        = Definitions.initguess;
maxIterations        = Definitions.maxIterations;
convergenceTolerance = Definitions.convergenceTolerance;
interval             = Definitions.interval;

% Set the initial conditions
phi0    = 0;
dphi0   = 0;
dphiInf = 1;

% Set ode solver options
opts = odeset('MaxStep', 1e-3);

% Generate the initial guesses for the shooting method
initializer = NaN(maxIterations+2, 1);
lam         = initializer;
dphi_inf    = initializer;
err         = initializer;
lam(1)      = lambdaGuesses(1);
lam(2)      = lambdaGuesses(2);
for n = 1:2
    x0          = [phi0 dphi0 lam(n)];
    [~, x]      = ode45(@(eta,x) hiemenzDoteq(eta,x), interval, x0, opts);
    dphi        = x(:,2);
    dphi_inf(n) = dphi(end);
    err(n)      = abs(dphi_inf(n)-dphiInf);
end

% Execute shoothing method to find the correct solution
for n = 3:maxIterations
    lam(n)      = abs(lam(n-1) + (lam(n-1)-lam(n-2))/(dphi_inf(n-1)-dphi_inf(n-2))*(1-dphi_inf(n-1)));
    x0          = [phi0 dphi0 lam(n)];
    [eta, x]    = ode45(@(eta,x) hiemenzDoteq(eta,x), interval, x0, opts);
    dphi        = x(:,2);
    dphi_inf(n) = dphi(end);
    err(n)    = abs(dphi_inf(n)-1);
    if err(n) < convergenceTolerance, break, end
end
phi   = x(:,1);
dphi  = x(:,2);
ddphi = x(:,3);

Base_Flow.eta   = eta;
Base_Flow.phi   = phi;
Base_Flow.dphi  = dphi;
Base_Flow.ddphi = ddphi;

end


function dx = hiemenzDoteq(~,x)

dx = zeros(numel(x),1);

phi   = x(1);
dphi  = x(2);
ddphi = x(3);

dx(1) = dphi;
dx(2) = ddphi;
dx(3) = dphi^2-phi*ddphi-1;

end