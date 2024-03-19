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

tic

%%

Tol = 1e-3;
EvalLim = 6e4;
p2c = 2;

f = @(z,p) D(z, p, p2c, 10);

y = [-0.5 0.5];
p = [1e-2 10];

SolNodes = [];

%%

n = 1;

x = [n n + 1];

sol = GAMES(f, [x; y; p], Tol, EvalLim);
SolNodes = [
    SolNodes;
    sol.ApproxNodes
    ];

%%

n = 2;

x = [n n + 1];

sol = GAMES(f, [x; y; p], Tol, EvalLim);
SolNodes = [
    SolNodes;
    sol.ApproxNodes
    ];

%%

n = 3;

x = [n n + 1];

sol = GAMES(f, [x; y; p], Tol, EvalLim);
SolNodes = [
    SolNodes;
    sol.ApproxNodes
    ];

%%

n = 4;

x = [n n + 1];

sol = GAMES(f, [x; y; p], Tol, EvalLim);
SolNodes = [
    SolNodes;
    sol.ApproxNodes
    ];

%%

n = 5;

x = [n n + 1];

sol = GAMES(f, [x; y; p], Tol, EvalLim);
SolNodes = [
    SolNodes;
    sol.ApproxNodes
    ];

%%

fprintf('Time elpased: %f s \n', toc);

ax = subplot(2,1,1);

hold on

for SolInd = [-1 1]
    SolPoints = SolNodes(SolNodes(:,end) == SolInd,1:end-1);

    scatter3( ...
        real(SolPoints(:,1)), ...
        imag(SolPoints(:,1)), ...
        SolPoints(:,2), ...
        5, ColorSet(SolInd+2,:), 'filled');
end

hold off

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$\omega / \omega_c$', FontSize=Axes_FontSize);
% ylabel(ax, '$y$', FontSize=Axes_FontSize);
zlabel(ax, '$k v / \omega_c$', FontSize=Axes_FontSize);
title(ax, '(a)', FontSize=Title_FontSize);

xlim([1 n + 1]);
ylim(y);
zlim(p);

view(ax,0,0);

tic

%%

Tol = 1e-3;
EvalLim = 6e4;
p2c = 4;

f = @(z,p) D(z, p, p2c, 10);

y = [-0.5 0.5];
p = [1e-2 10];

SolNodes = [];

%%

n = 1;

x = [n n + 1];

sol = GAMES(f, [x; y; p], Tol, EvalLim);
SolNodes = [
    SolNodes;
    sol.ApproxNodes
    ];

%%

n = 2;

x = [n n + 1];

sol = GAMES(f, [x; y; p], Tol, EvalLim);
SolNodes = [
    SolNodes;
    sol.ApproxNodes
    ];

%%

n = 3;

x = [n n + 1];

sol = GAMES(f, [x; y; p], Tol, EvalLim);
SolNodes = [
    SolNodes;
    sol.ApproxNodes
    ];

%%

n = 4;

x = [n n + 1];

sol = GAMES(f, [x; y; p], Tol, EvalLim);
SolNodes = [
    SolNodes;
    sol.ApproxNodes
    ];

%%

n = 5;

x = [n n + 1];

sol = GAMES(f, [x; y; p], Tol, EvalLim);
SolNodes = [
    SolNodes;
    sol.ApproxNodes
    ];

%%

fprintf('Time elpased: %f s \n', toc);

ax = subplot(2,1,2);

hold on

for SolInd = [-1 1]
    SolPoints = SolNodes(SolNodes(:,end) == SolInd,1:end-1);

    scatter3( ...
        real(SolPoints(:,1)), ...
        imag(SolPoints(:,1)), ...
        SolPoints(:,2), ...
        5, ColorSet(SolInd+2,:), 'filled');
end

hold off

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$\omega / \omega_c$', FontSize=Axes_FontSize);
% ylabel(ax, '$y$', FontSize=Axes_FontSize);
zlabel(ax, '$k v / \omega_c$', FontSize=Axes_FontSize);
title(ax, '(b)', FontSize=Title_FontSize);

xlim([1 n + 1]);
ylim(y);
zlim(p);

view(ax,0,0);

function value = D(omega, k, p2c, n_max)
    value = 1 - besseli(0, k.^2) .* exp(-k.^2);

    for n=1:n_max
        value = value - (besseli(n, k.^2) ./ (omega - n) + besseli(-n, k.^2) ./ (omega + n)) .* omega .* exp(-k.^2);
    end

    value = 1 + p2c^2 ./ k.^2 .* value;
end