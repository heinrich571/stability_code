function Base_Flow = spectralGetBaseFlow(Definitions)

% Get data from user input / fixed user input
lambdaGuesses        = Definitions.initguess;
maxIterations        = Definitions.maxIterations;
convergenceTolerance = Definitions.convergenceTolerance;
interval             = flip(Definitions.interval);

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

end

% === Supporting Functions === %
% Solving the differential equation
function 
