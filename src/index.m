%% ENAE 641 - Homework 3
% Vai Srivastava

% clear history
clear;
clc;
close all;

% resolve paths from this file, so the script runs from any working folder
root = fileparts(fileparts(mfilename("fullpath")));
textDir = fullfile(root, "outputs", "text");
figDir = fullfile(root, "outputs", "figures");
if ~isfolder(textDir), mkdir(textDir); end
if ~isfolder(figDir), mkdir(figDir); end

% load provided data
load(fullfile(root, "data", "F16.mat"), "A", "B", "C");
n = size(A, 1);
stateNames = ["speed"; "angle of attack"; "pitch angle"; "pitch rate"];
sys = ss(A, B, C, 0);

% categorical series colours (fixed order) and chart ink
colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];
muted = "#52514e";

%% Problem 1a: real block modal form
[V, D] = eig(A);
lambda = diag(D);

% keep the +j*omega member of each conjugate pair, ordered slow -> fast
% (phugoid, then short period)
keep = imag(lambda) > 0;
lambdaModes = lambda(keep);
vModes = V(:, keep);
[~, order] = sort(abs(lambdaModes));
lambdaModes = lambdaModes(order);
vModes = vModes(:, order);

% scale each eigenvector so that its largest entry is exactly 1
for k = 1:2
    [~, iMax] = max(abs(vModes(:, k)));
    vModes(:, k) = vModes(:, k) / vModes(iMax, k);
end
vPh = vModes(:, 1);
vSp = vModes(:, 2);

% x = Pm z, with Pm = [Re(v_ph), Im(v_ph), Re(v_sp), Im(v_sp)]
Pm = [real(vPh), imag(vPh), real(vSp), imag(vSp)];
Ap = Pm \ A * Pm;
Ap(abs(Ap) < 1e-12) = 0;
Bp = Pm \ B;

sigma = real(lambdaModes);
omega = imag(lambdaModes);
wn = abs(lambdaModes);
zeta = -sigma ./ wn;
period = 2 * pi ./ omega;
modes = table(["phugoid"; "short period"], sigma, omega, wn, zeta, period, ...
    VariableNames=["mode", "sigma", "omega", "omega_n", "zeta", "period_s"]);

report(textDir, "s01a", join([
    "eigenvalues of A:"
    show(lambda)
    ""
    "modes (lambda = sigma +/- j omega):"
    show(modes)
    ""
    "max Re(lambda) = " + num2str(max(real(lambda)), 6) + " < 0  ->  asymptotically stable"
    ""
    "Pm = [Re(v_ph), Im(v_ph), Re(v_sp), Im(v_sp)]:"
    show(Pm)
    ""
    "A_p = Pm^-1 A Pm:"
    show(Ap)
], newline));

%% Problem 1b: single-mode initial conditions
polarPh = table(stateNames, abs(vPh), rad2deg(angle(vPh)), VariableNames=["state", "r", "phi_deg"]);
polarSp = table(stateNames, abs(vSp), rad2deg(angle(vSp)), VariableNames=["state", "r", "phi_deg"]);

% x0 in span{Re(v), Im(v)} excites only that mode; take x0 = Re(v)
x0Ph = real(vPh);
x0Sp = real(vSp);

tPh = linspace(0, 400, 2001)';
tSp = linspace(0, 6, 1201)';
[~, tPh, xPh] = initial(sys, x0Ph, tPh);
[~, tSp, xSp] = initial(sys, x0Sp, tSp);
zPh = (Pm \ xPh')';
zSp = (Pm \ xSp')';

report(textDir, "s01b", join([
    "phugoid eigenvector v_ph (polar form, relative to speed):"
    show(polarPh)
    ""
    "short-period eigenvector v_sp (polar form, relative to pitch rate):"
    show(polarSp)
    ""
    "phugoid x0 = Re(v_ph):"
    show(x0Ph')
    "z0 = Pm^-1 x0:"
    show((Pm \ x0Ph)')
    "max |z3|, |z4| over the phugoid response:  " + num2str(max(abs(zPh(:, 3:4)), [], "all"), 3)
    ""
    "short-period x0 = Re(v_sp):"
    show(x0Sp')
    "z0 = Pm^-1 x0:"
    show((Pm \ x0Sp)')
    "max |z1|, |z2| over the short-period response:  " + num2str(max(abs(zSp(:, 1:2)), [], "all"), 3)
], newline));

fig = figure(Visible="off", Units="centimeters", Position=[0, 0, 18, 12]);
layout = tiledlayout(fig, 2, 2, TileSpacing="compact", Padding="compact");
cases = {tPh, xPh, zPh, "phugoid: x_0 = Re\{v_{ph}\}"; tSp, xSp, zSp, "short period: x_0 = Re\{v_{sp}\}"};
for c = 1:2
    [t, x, z, label] = cases{c, :};

    ax = nexttile(layout, c);
    hold(ax, "on");
    for i = 1:n
        plot(ax, t, x(:, i), Color=colors(i), LineWidth=1.5);
    end
    styleAxes(ax, muted);
    title(ax, label);
    ylabel(ax, "state x_i(t)");
    if c == 1
        legend(ax, stateNames, Location="northeast", Box="off");
    end

    ax = nexttile(layout, c + 2);
    hold(ax, "on");
    for i = 1:n
        plot(ax, t, z(:, i), Color=colors(i), LineWidth=1.5);
    end
    styleAxes(ax, muted);
    xlabel(ax, "time t (s)");
    ylabel(ax, "modal coordinate z_i(t)");
    if c == 1
        legend(ax, "z_" + (1:n), Location="northeast", Box="off");
    end
end
saveFigure(fig, figDir, "s01b");

%% Problem 1c: modal B matrix
% rows: phugoid block (z1, z2), short-period block (z3, z4); columns: throttle, elevator
blocks = {1:2, 3:4};
blockNorm = zeros(2, 2);
projNorm = zeros(2, 2);
for k = 1:2
    for j = 1:2
        blockNorm(k, j) = norm(Bp(blocks{k}, j));
        % part of b_j lying in mode k's plane (independent of eigenvector scaling)
        projNorm(k, j) = norm(Pm(:, blocks{k}) * Bp(blocks{k}, j));
    end
end
inputs = ["throttle", "elevator"];
blockTable = array2table(blockNorm, VariableNames=inputs, RowNames=["phugoid", "short period"]);
projTable = array2table(projNorm, VariableNames=inputs, RowNames=["phugoid", "short period"]);

report(textDir, "s01c", join([
    "B_p = Pm^-1 B (columns: throttle, elevator):"
    show(Bp)
    ""
    "modal-block norms ||B_p(block, j)||:"
    show(blockTable)
    "ratios (larger / smaller): throttle " + num2str(blockNorm(1, 1) / blockNorm(2, 1), 3) ...
        + ", elevator " + num2str(blockNorm(2, 2) / blockNorm(1, 2), 3)
    ""
    "projection of each input direction onto each mode's plane, ||Pm_k B_p,k(:, j)||:"
    show(projTable)
    "ratios (larger / smaller): throttle " + num2str(projNorm(1, 1) / projNorm(2, 1), 3) ...
        + ", elevator " + num2str(projNorm(2, 2) / projNorm(1, 2), 3)
], newline));

%% Problem 4a: Lyapunov matrix
% lyap(A, Q) solves A*X + X*A' + Q = 0, so pass A' to solve A'*P + P*A + Q = 0
Q = eye(n);
P = lyap(A', Q);
eigP = eig(P);
[~, cholFlag] = chol(P);

report(textDir, "s04a", join([
    "P = lyap(A', I):"
    show(P)
    ""
    "eigenvalues of P:"
    show(eigP')
    "min eigenvalue = " + num2str(min(eigP), 6) + " > 0,  chol flag = " + cholFlag ...
        + "  ->  P is positive definite"
], newline));

%% Problem 4b: exponential convergence bound
rho = max(eigP);
eta = min(eigP);
nu = min(eig(Q));
gamma = sqrt(rho / eta);
lambdaBound = nu / (2 * rho);
slowestRate = min(-real(eig(A)));

tB = (0:0.01:600)';
gainA = arrayfun(@(t) norm(expm(A * t)), tB);
bound = gamma * exp(-lambdaBound * tB);
[peakGain, iPeak] = max(gainA);

% where the conservatism comes from
phPlane = orth(Pm(:, 1:2));
rhoPh = max(eig(phPlane' * P * phPlane));
reV = Pm(:, 1);
imV = Pm(:, 2);
phAngle = acosd(reV' * imV / (norm(reV) * norm(imV)));
minM = max(gainA .* exp(lambdaBound * tB));

report(textDir, "s04b", join([
    "rho = lambda_max(P) = " + num2str(rho, 6)
    "eta = lambda_min(P) = " + num2str(eta, 6)
    "nu  = lambda_min(Q) = " + num2str(nu, 6)
    "gamma  = sqrt(rho / eta) = " + num2str(gamma, 6)
    "lambda = nu / (2 rho)   = " + num2str(lambdaBound, 6) + " 1/s"
    ""
    "slowest decay rate of A, min |Re(lambda)| = " + num2str(slowestRate, 6) + " 1/s"
    "max_t ||e^(At)||_2 = " + num2str(peakGain, 6) + " at t = " + num2str(tB(iPeak), 4) + " s"
    "bound >= ||e^(At)||_2 on the whole grid: " + string(all(bound >= gainA))
    ""
    "lambda_max of P restricted to the phugoid plane = " + num2str(rhoPh, 4) ...
        + "  (1 / (2 |sigma_1|) = " + num2str(1 / (2 * abs(sigma(1))), 4) + ")"
    "|Re(v_ph)| = " + num2str(norm(reV), 3) + ", |Im(v_ph)| = " + num2str(norm(imV), 3) ...
        + ", angle between them = " + num2str(phAngle, 3) + " deg"
    "smallest M with ||e^(At)||_2 <= M e^(-lambda t): " + num2str(minM, 4)
    "|sigma_3| / |sigma_1| = " + num2str(sigma(2) / sigma(1), 4)
    "modal bound: ||x(t)|| <= cond(Pm) e^(sigma_1 t) ||x0||, cond(Pm) = " + num2str(cond(Pm), 4)
], newline));

fig = figure(Visible="off", Units="centimeters", Position=[0, 0, 16, 8]);
ax = axes(fig);
hold(ax, "on");
plot(ax, tB, gainA, Color=colors(1), LineWidth=1.5);
plot(ax, tB, bound, Color=colors(2), LineWidth=1.5);
ax.YScale = "log";
styleAxes(ax, muted);
xlabel(ax, "time t (s)");
ylabel(ax, "||x(t)||_2 / ||x_0||_2");
legend(ax, ["worst case ||e^{At}||_2", "Lyapunov bound \gamma e^{-\lambda t}"], Location="northeast", Box="off");
saveFigure(fig, figDir, "s04b");

%% Problem 5a: eigenvalues of A + A'
M = A + A';
eigM = eig(M);

report(textDir, "s05a", join([
    "eigenvalues of A + A':"
    show(eigM')
    "number of positive eigenvalues: " + sum(eigM > 0)
], newline));

%% Problem 5b: ||x||^2 vs x'Px along solutions
% x0 = w: ||x||^2 grows at t = 0 (largest eigenvalue of A + A')
[W, E] = eig(M);
[~, iTop] = max(diag(E));
w = W(:, iTop);
growth0 = w' * M * w;

% x0 = xw: worst case at the peak of ||e^(At)||_2 (top right singular vector)
[~, ~, Vsvd] = svd(expm(A * tB(iPeak)));
xw = Vsvd(:, 1);

t5 = (0:0.005:60)';
ics = [w, xw];
normRatio = zeros(numel(t5), 2);
vRatio = zeros(numel(t5), 2);
for c = 1:2
    [~, ~, x5] = initial(sys, ics(:, c), t5);
    normRatio(:, c) = vecnorm(x5, 2, 2) / norm(ics(:, c));
    vRatio(:, c) = sum((x5 * P) .* x5, 2) / (ics(:, c)' * P * ics(:, c));
end
[peakNorm, iPk] = max(normRatio);
iBack = find(normRatio(:, 1) < 1 & t5 > t5(iPk(1)), 1);

report(textDir, "s05b", join([
    "x0 = w, the unit eigenvector of A + A' for its largest eigenvalue:"
    show(w')
    "d/dt ||x||^2 at t = 0: w'(A + A')w = " + num2str(growth0, 6) + " > 0"
    "max ||x(t)|| / ||x0|| = " + num2str(peakNorm(1), 6) + " at t = " + num2str(t5(iPk(1)), 4) ...
        + " s (back below 1 at t = " + num2str(t5(iBack), 4) + " s)"
    ""
    "x0 = x_w, the unit vector maximizing ||e^(At)x0|| at t = " + num2str(tB(iPeak), 4) + " s:"
    show(xw')
    "max ||x(t)|| / ||x0|| = " + num2str(peakNorm(2), 6) + " at t = " + num2str(t5(iPk(2)), 4) + " s"
    ""
    "V = x'Px nonincreasing along both solutions: " + string(all(diff(vRatio) <= 0, "all"))
], newline));

fig = figure(Visible="off", Units="centimeters", Position=[0, 0, 18, 7.5]);
layout = tiledlayout(fig, 1, 3, TileSpacing="compact", Padding="compact");
% norm over 60 s, the same norm over the first second, and V over 60 s
panels = {normRatio, 60, "||x(t)||_2 / ||x_0||_2"; normRatio, 1, "||x(t)||_2 / ||x_0||_2, first 1 s"; ...
    vRatio, 60, "V(t) / V(0),  V = x^T P x"};
for c = 1:3
    ax = nexttile(layout);
    hold(ax, "on");
    yline(ax, 1, Color="#c3c2b7", LineWidth=0.8, HandleVisibility="off");
    inWindow = t5 <= panels{c, 2};
    for i = 1:2
        plot(ax, t5(inWindow), panels{c, 1}(inWindow, i), Color=colors(i), LineWidth=1.5);
    end
    styleAxes(ax, muted);
    xlabel(ax, "time t (s)");
    title(ax, panels{c, 3});
    if c == 1
        legend(ax, ["x_0 = w", "x_0 = x_w"], Location="southwest", Box="off");
    end
end
saveFigure(fig, figDir, "s05b");

%% Helper functions
function out = show(value)
    % display text for a value, without trailing blank lines
    out = strip(formattedDisplayText(value, NumericFormat="shortg", ...
        LineSpacing="compact", SuppressMarkup=true), "right");
end

function report(dir, name, txt)
    % print a result and save it for the report
    fprintf("\n------\n %s\n------\n%s\n", name, txt);
    writelines(txt, fullfile(dir, name + ".txt"));
end

function styleAxes(ax, ink)
    % recessive axes and grid
    ax.Box = "off";
    ax.FontSize = 9;
    ax.TickDir = "out";
    ax.LineWidth = 0.6;
    ax.XColor = ink;
    ax.YColor = ink;
    ax.XGrid = "on";
    ax.YGrid = "on";
    ax.GridColor = "#c3c2b7";
    ax.GridAlpha = 0.5;
end

function saveFigure(fig, dir, name)
    % export a figure for the report
    exportgraphics(fig, fullfile(dir, name + ".png"), Resolution=200);
    close(fig);
end
