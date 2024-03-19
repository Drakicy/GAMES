set(0, defaulttextinterpreter='latex')

f = @(z,p) z - p + exp(10 * (z + p));

x = [-2 2];
y = x;
p = x;

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

Tol = 1e-3;

%%

EvalLim = 5e3;

ax = subplot(2,2,1);

tic

sol = GAMES(f, [x; y; p], Tol, EvalLim);

fprintf('Time elpased: %f s \n', toc);

hold on

for SolInd = [-1 1 0]
    SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == SolInd,1:end-1);

    scatter3( ...
        real(SolPoints(:,1)), ...
        imag(SolPoints(:,1)), ...
        SolPoints(:,2), ...
        5, ColorSet(SolInd+2,:), 'filled');
end

hold off

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$x$', FontSize=Axes_FontSize);
ylabel(ax, '$y$', FontSize=Axes_FontSize);
zlabel(ax, '$p$', FontSize=Axes_FontSize);
title(ax, '(a)', FontSize=Title_FontSize);

xlim(x);
ylim(y);
zlim(p);

view(ax,45,30);

%%

EvalLim = 1e4;

ax = subplot(2,2,2);

tic

sol = GAMES(f, [x; y; p], Tol, EvalLim);

fprintf('Time elpased: %f s \n', toc);

hold on

for SolInd = [-1 1 0]
    SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == SolInd,1:end-1);

    scatter3( ...
        real(SolPoints(:,1)), ...
        imag(SolPoints(:,1)), ...
        SolPoints(:,2), ...
        5, ColorSet(SolInd+2,:), 'filled');
end

hold off

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$x$', FontSize=Axes_FontSize);
ylabel(ax, '$y$', FontSize=Axes_FontSize);
zlabel(ax, '$p$', FontSize=Axes_FontSize);
title(ax, '(b)', FontSize=Title_FontSize);

xlim(x);
ylim(y);
zlim(p);

view(ax,45,30);

%%

EvalLim = 2e4;

ax = subplot(2,2,3);

tic

sol = GAMES(f, [x; y; p], Tol, EvalLim);

fprintf('Time elpased: %f s \n', toc);

hold on

for SolInd = [-1 1 0]
    SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == SolInd,1:end-1);

    scatter3( ...
        real(SolPoints(:,1)), ...
        imag(SolPoints(:,1)), ...
        SolPoints(:,2), ...
        5, ColorSet(SolInd+2,:), 'filled');
end

hold off

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$x$', FontSize=Axes_FontSize);
ylabel(ax, '$y$', FontSize=Axes_FontSize);
zlabel(ax, '$p$', FontSize=Axes_FontSize);
title(ax, '(c)', FontSize=Title_FontSize);

xlim(x);
ylim(y);
zlim(p);

view(ax,45,30);

%%

EvalLim = 4e4;

ax = subplot(2,2,4);

tic

sol = GAMES(f, [x; y; p], Tol, EvalLim);

fprintf('Time elpased: %f s \n', toc);

hold on

for SolInd = [-1 1 0]
    SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == SolInd,1:end-1);

    scatter3( ...
        real(SolPoints(:,1)), ...
        imag(SolPoints(:,1)), ...
        SolPoints(:,2), ...
        5, ColorSet(SolInd+2,:), 'filled');
end

hold off

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$x$', FontSize=Axes_FontSize);
ylabel(ax, '$y$', FontSize=Axes_FontSize);
zlabel(ax, '$p$', FontSize=Axes_FontSize);
title(ax, '(d)', FontSize=Title_FontSize);

xlim(x);
ylim(y);
zlim(p);

view(ax,45,30);