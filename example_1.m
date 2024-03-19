set(0, defaulttextinterpreter='latex')

f = @(z,a) (z-a) ./ (z + a);

x = [-2 2];
y = x;

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

a = 0.9;
Tol = 1e-3;

%%

EvalLim = 5e1;

ax = subplot(2,2,1);

tic

sol = GAMES(@(x) f(x,a), [x; y], Tol, EvalLim);

fprintf('Time elpased: %f s \n', toc);

ResizedDT = sol.DT;
ResizedDT.Points = ResizedDT.Points .* sol.GridNorm(:,2).' + sol.GridNorm(:,1).';

hold on

triplot(ResizedDT, '--k');

for SolInd = [-1 1]
    SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == SolInd,1:end-1);

    scatter(ax, ...
        real(SolPoints), ...
        imag(SolPoints), ...
        25, ColorSet(SolInd+2,:), 'filled');
end

hold off

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$x$', FontSize=Axes_FontSize);
ylabel(ax, '$y$', FontSize=Axes_FontSize);
title(ax, '(a)', FontSize=Title_FontSize);

xlim(x);
ylim(y);

%%

EvalLim = 1e2;

ax = subplot(2,2,2);

tic

sol = GAMES(@(x) f(x,a), [x; y], Tol, EvalLim);

fprintf('Time elpased: %f s \n', toc);

ResizedDT = sol.DT;
ResizedDT.Points = ResizedDT.Points .* sol.GridNorm(:,2).' + sol.GridNorm(:,1).';

hold on

triplot(ResizedDT, '--k');

for SolInd = [-1 1 0]
    SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == SolInd,1:end-1);

    scatter(ax, ...
        real(SolPoints), ...
        imag(SolPoints), ...
        25, ColorSet(SolInd+2,:), 'filled');
end

hold off

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$x$', FontSize=Axes_FontSize);
ylabel(ax, '$y$', FontSize=Axes_FontSize);
title(ax, '(b)', FontSize=Title_FontSize);

xlim(x);
ylim(y);

%%

EvalLim = 2e2;

ax = subplot(2,2,3);

tic

sol = GAMES(@(x) f(x,a), [x; y], Tol, EvalLim);

fprintf('Time elpased: %f s \n', toc);

ResizedDT = sol.DT;
ResizedDT.Points = ResizedDT.Points .* sol.GridNorm(:,2).' + sol.GridNorm(:,1).';

hold on

triplot(ResizedDT, '--k');

for SolInd = [-1 1 0]
    SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == SolInd,1:end-1);

    scatter(ax, ...
        real(SolPoints), ...
        imag(SolPoints), ...
        25, ColorSet(SolInd+2,:), 'filled');
end

hold off

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$x$', FontSize=Axes_FontSize);
ylabel(ax, '$y$', FontSize=Axes_FontSize);
title(ax, '(c)', FontSize=Title_FontSize);

xlim(x);
ylim(y);

%%

EvalLim = 5e2;

ax = subplot(2,2,4);

tic

sol = GAMES(@(x) f(x,a), [x; y], Tol, EvalLim);

fprintf('Time elpased: %f s \n', toc);

ResizedDT = sol.DT;
ResizedDT.Points = ResizedDT.Points .* sol.GridNorm(:,2).' + sol.GridNorm(:,1).';

hold on

triplot(ResizedDT, '--k');

for SolInd = [-1 1 0]
    SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == SolInd,1:end-1);

    scatter(ax, ...
        real(SolPoints), ...
        imag(SolPoints), ...
        25, ColorSet(SolInd+2,:), 'filled');
end

hold off

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$x$', FontSize=Axes_FontSize);
ylabel(ax, '$y$', FontSize=Axes_FontSize);
title(ax, '(d)', FontSize=Title_FontSize);

xlim(x);
ylim(y);

%%

fig = figure(Position=[0 0 1000 500]);

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
EvalLim = 5e2;

%%

a = 1e-1;

ax = subplot(1,2,1);

tic

sol = GAMES(@(x) f(x,a), [x; y], Tol, EvalLim);

fprintf('Time elpased: %f s \n', toc);

ResizedDT = sol.DT;
ResizedDT.Points = ResizedDT.Points .* sol.GridNorm(:,2).' + sol.GridNorm(:,1).';

hold on

triplot(ResizedDT, '--k');

for SolInd = [-1 1 0]
    SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == SolInd,1:end-1);

    scatter(ax, ...
        real(SolPoints), ...
        imag(SolPoints), ...
        25, ColorSet(SolInd+2,:), 'filled');
end

hold off

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$x$', FontSize=Axes_FontSize);
ylabel(ax, '$y$', FontSize=Axes_FontSize);
title(ax, '(a)', FontSize=Title_FontSize);

xlim(x);
ylim(y);

%%

a = 1e-2;

ax = subplot(1,2,2);

tic

sol = GAMES(@(x) f(x,a), [x; y], Tol, EvalLim);

fprintf('Time elpased: %f s \n', toc);

ResizedDT = sol.DT;
ResizedDT.Points = ResizedDT.Points .* sol.GridNorm(:,2).' + sol.GridNorm(:,1).';

hold on

triplot(ResizedDT, '--k');

for SolInd = [-1 1 0]
    SolPoints = sol.ApproxNodes(sol.ApproxNodes(:,end) == SolInd,1:end-1);

    scatter(ax, ...
        real(SolPoints), ...
        imag(SolPoints), ...
        25, ColorSet(SolInd+2,:), 'filled');
end

hold off

set(ax, FontSize=Global_FontSize);
xlabel(ax, '$x$', FontSize=Axes_FontSize);
ylabel(ax, '$y$', FontSize=Axes_FontSize);
title(ax, '(b)', FontSize=Title_FontSize);

xlim(x);
ylim(y);