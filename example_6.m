set(0, defaulttextinterpreter='latex')

fig = figure(Position=[0 0 1000 1000]);

Global_FontSize = 18;
Axes_FontSize = 25;
Title_FontSize = 22;

ColorSet = ...
    [
    0.8500 0.3250 0.0980;
    0 0 0;
    0 0.4470 0.7410
    ];

%%

Tol = 1e-3;
EvalLim = 1e6;

m = [1 1 1836];
n = [1 1 2];
T = [1 1 1];

f = @(z,p) D(z .* exp(1i * pi / 4), p, m, n, T);

x = [1e-1 2];
y = [-5 0];
p = [0 2];

tic

sol = GAMES(f, [x; y; p], Tol, EvalLim);

fprintf('Time elpased: %f s \n', toc);

%%

ax = subplot(2,2,1);

SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == 1,1:end-1);
SolPoints(:,1) = SolPoints(:,1) .* exp(1i * pi / 4) / sqrt(2);

scatter3( ...
    real(SolPoints(:,1)), ...
    imag(SolPoints(:,1)), ...
    SolPoints(:,2), ...
    5, ColorSet(3,:), 'filled');

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$\Omega/\omega_{pe}$', FontSize=Axes_FontSize);
% ylabel(ax, '$y$', FontSize=Axes_FontSize);
zlabel(ax, '$k/k_{Dc}$', FontSize=Axes_FontSize);
title(ax, '(a)', FontSize=Title_FontSize);

view(ax,0,0);

xlim([0 3.5 / sqrt(2)]);
% ylim([-3 / sqrt(2) 0]);
zlim([0 2]);

%%

ax = subplot(2,2,2);

SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == 1,1:end-1);
SolPoints(:,1) = SolPoints(:,1) .* exp(1i * pi / 4) / sqrt(2);

scatter3( ...
    real(SolPoints(:,1)), ...
    imag(SolPoints(:,1)), ...
    SolPoints(:,2), ...
    5, ColorSet(3,:), 'filled');

set(ax, FontSize=Global_FontSize);
set(ax, 'YDir', 'reverse');
% xlabel(ax, '$x$', FontSize=Axes_FontSize);
ylabel(ax, '$\gamma/\omega_{pe}$', FontSize=Axes_FontSize);
zlabel(ax, '$k/k_{Dc}$', FontSize=Axes_FontSize);
title(ax, '(b)', FontSize=Title_FontSize);

view(ax,90,0);

% xlim([0 3.5 / sqrt(2)]);
ylim([-3 / sqrt(2) 0]);
zlim([0 2]);

%%

Tol = 1e-3;
EvalLim = 1e6;

m = [1 1 1836];
n = [1 1 2];
T = [1 10 1];

f = @(z,p) D(z .* exp(1i * pi / 4), p, m, n, T);

x = [1e-1 2];
y = [-5 0];
p = [0 2];

tic

sol = GAMES(f, [x; y; p], Tol, EvalLim);

fprintf('Time elpased: %f s \n', toc);

%%

ax = subplot(2,2,3);

SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == 1,1:end-1);
SolPoints(:,1) = SolPoints(:,1) .* exp(1i * pi / 4) / sqrt(2);

scatter3( ...
    real(SolPoints(:,1)), ...
    imag(SolPoints(:,1)), ...
    SolPoints(:,2), ...
    5, ColorSet(3,:), 'filled');

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$\Omega/\omega_{pe}$', FontSize=Axes_FontSize);
% ylabel(ax, '$y$', FontSize=Axes_FontSize);
zlabel(ax, '$k/k_{Dc}$', FontSize=Axes_FontSize);
title(ax, '(c)', FontSize=Title_FontSize);

view(ax,0,0);

xlim([0 3.5 / sqrt(2)]);
% ylim([-3 / sqrt(2) 0]);
zlim([0 2]);

%%

ax = subplot(2,2,4);

SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == 1,1:end-1);
SolPoints(:,1) = SolPoints(:,1) .* exp(1i * pi / 4) / sqrt(2);

scatter3( ...
    real(SolPoints(:,1)), ...
    imag(SolPoints(:,1)), ...
    SolPoints(:,2), ...
    5, ColorSet(3,:), 'filled');

set(ax, FontSize=Global_FontSize);
set(ax, 'YDir', 'reverse');
xlabel(ax, '$x$', FontSize=Axes_FontSize);
ylabel(ax, '$\gamma/\omega_{pe}$', FontSize=Axes_FontSize);
% zlabel(ax, '$k/k_{Dc}$', FontSize=Axes_FontSize);
title(ax, '(d)', FontSize=Title_FontSize);

view(ax,90,0);

% xlim([0 3.5 / sqrt(2)]);
ylim([-3 / sqrt(2) 0]);
zlim([0 2]);

function value = D(omega, k, m, n, T)
    m = reshape(m,1,1,[]);
    n = reshape(n,1,1,[]);
    T = reshape(T,1,1,[]);
    omegaInd = (omega == 0);
    kInd = (k == 0);

    value = zeros(size(omega));
    value(~omegaInd & kInd) = 1 - sum(n ./ m,3) ./ omega(~omegaInd & kInd).^2;
    value(~kInd) = 1 - 1 ./ (2 * k(~kInd).^2) .* sum(n ./ T .* pdfZDiv(omega(~kInd) ./ k(~kInd) .* sqrt(m ./ (2 * T))),3);
end