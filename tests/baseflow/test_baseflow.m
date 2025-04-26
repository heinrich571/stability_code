%% Fresh Start

close all
clear all
clc

path_manager('add')
% addpath('../../baseflow_definition/')
% addpath('../../grid_generation/')


%% Define Parameters

ymax_vec = 300; % values of maximum practical farfield distance from wall


%% Load Reference Solution Data

refsol = readtable("./hiemenz_flow_data.csv");


%% Base Flow Calculation

Base_Flow_Definitions.initguess            = 1.23258765682022 + [-1 1]*1e-5;
Base_Flow_Definitions.maxIterations        = 1e2;
Base_Flow_Definitions.convergenceTolerance = 1e-6;

for i = 1 : length(ymax_vec)
    ymax = ymax_vec(i);

    Base_Flow_Definitions.interval = flip([refsol.eta ; ymax]);

    Base_Flow(i) = get_base_flow(Base_Flow_Definitions);

    Base_Flow(i).eta   = Base_Flow(i).eta(1:end-1);
    Base_Flow(i).phi   = Base_Flow(i).phi(1:end-1);
    Base_Flow(i).dphi  = Base_Flow(i).dphi(1:end-1);
    Base_Flow(i).ddphi = Base_Flow(i).ddphi(1:end-1);
end


%% Errors Calculation

for i = 1 : length(Base_Flow)
    phi   = Base_Flow(i).phi;
    dphi  = Base_Flow(i).dphi;
    ddphi = Base_Flow(i).ddphi;

    err(i).eta   = Base_Flow(i).eta;
    err(i).phi   = abs(phi   - refsol.phi        );
    err(i).dphi  = abs(dphi  - refsol.dphi_deta  );
    err(i).ddphi = abs(ddphi - refsol.d2phi_deta2);
end



%% Plots

howarth_label = "Howarth, 1939";

% phi(eta)
figure("Name", "phi(eta)")
grid on
hold on
plot(refsol.eta, refsol.phi, "dk", "MarkerSize", 10, "DisplayName", howarth_label)
for i = 1 : length(Base_Flow)
    ymax = ymax_vec(i);

    plot(Base_Flow(i).eta, Base_Flow(i).phi, "-xb", "DisplayName", strcat("Present work, $\eta_{max} = ", num2str(ymax), "$"))
end
legend("Interpreter", "latex", "Location", "northwest")
xlabel("$\eta$", "Interpreter", "latex")
ylabel("$\phi(\eta)$", "Interpreter", "latex")
ax = gca;
ax.YColor = 'b';
yyaxis right
ax = gca;
ax.YColor = 'r';
for i = 1 : length(err)
    ymax = ymax_vec(i);

    plot(err(i).eta, log10(err(i).phi), ":xr", "HandleVisibility", "off")
end
ylabel("$\mathrm{log}_{10}\left(\epsilon_{\phi\left(\eta\right)}\right)$", "Interpreter", "latex")


% dphi(eta)
figure("Name", "dphi(eta)")
grid on
hold on
plot(refsol.eta, refsol.dphi_deta, "dk", "MarkerSize", 10, "DisplayName", howarth_label)
for i = 1 : length(Base_Flow)
    ymax = ymax_vec(i);

    plot(Base_Flow(i).eta, Base_Flow(i).dphi, "-xb", "DisplayName", strcat("Present work, $\eta_{max} = ", num2str(ymax), "$"))
end
legend("Interpreter", "latex", "Location", "northwest")
xlabel("$\phi'(\eta)$", "Interpreter", "latex")
ylabel("$\eta$", "Interpreter", "latex")
ax = gca;
ax.YColor = 'b';
ylim([0 1.2])
yyaxis right
ax = gca;
ax.YColor = 'r';
for i = 1 : length(err)
    ymax = ymax_vec(i);

    plot(err(i).eta, log10(err(i).dphi), ":xr", "HandleVisibility", "off")
end
ylabel("$\mathrm{log}_{10}\left(\epsilon_{\phi'\left(\eta\right)}\right)$", "Interpreter", "latex")


% ddphi(eta)
figure("Name", "ddphi(eta)")
grid on
hold on
plot(refsol.eta, refsol.d2phi_deta2, "dk", "MarkerSize", 10, "DisplayName", howarth_label)
for i = 1 : length(Base_Flow)
    ymax = ymax_vec(i);

    plot(Base_Flow(i).eta, Base_Flow(i).ddphi, "-xb", "DisplayName", strcat("Present work, $\eta_{max} = ", num2str(ymax), "$"))
end
legend("Interpreter", "latex", "Location", "northwest")
xlabel("$\phi''(\eta)$", "Interpreter", "latex")
ylabel("$\eta$", "Interpreter", "latex")
ax = gca;
ax.YColor = 'b';
yyaxis right
ax = gca;
ax.YColor = 'r';
for i = 1 : length(err)
    ymax = ymax_vec(i);

    plot(err(i).eta, log10(err(i).ddphi), ":xr", "HandleVisibility", "off")
end
ylabel("$\mathrm{log}_{10}\left(\epsilon_{\phi''\left(\eta\right)}\right)$", "Interpreter", "latex")



%% Cleanup Paths

path_manager('remove')
% rmpath('..\..\baseflow_definition\')
% rmpath('..\..\grid_generation\')

