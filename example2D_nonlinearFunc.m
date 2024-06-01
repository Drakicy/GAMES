clear
set(0, defaulttextinterpreter='latex')

Global_FontSize = 18;
Axes_FontSize = 25;
Title_FontSize = 22;

Title_Set = ["(a)" "(b)" "(c)" "(d)"];
ColorSet = ...
    [
    0.8500 0.3250 0.0980;
    0 0 0;
    0 0.4470 0.7410
    ];

%%

f = @(z) (exp(4 * pi * (z + (sqrt(3) + 1i * sqrt(2)) / 4)) - 1) ./ (exp(4 * pi * (z + (sqrt(3) + 1i * sqrt(2)) / 4)) + 1);

x = [-1 1];
y = x;

Tol = 1e-3;

%%

PointNumLim = [5e1 1e2 2e2 floor(1 / prod(Tol) / Tol(1))];

%%

fig = figure(Position=[0 0 1000 1000]);

t = 0;

for i = 1:4
    tic

    %%

    if exist('sol', 'var')
        sol.PointNumLim = PointNumLim(i);
        sol.fitTriang;
    else
        sol = GAMES(f, [x; y]', Tol, PointNumLim(i), Display=0);
    end

    %%

    t = t + toc;

    fprintf('Step <strong>#%i</strong> done!\nOverall time: <strong>%f sec</strong>\nPoints number: <strong>%i</strong>\n', i, t, sol.PointNum);

    %%

    ax = subplot(2,2,i);

    hold(ax, 'on')
    
    TR = triangulation(sol.DT(:,:), sol.GridNorm(1,:) + sol.GridNorm(2,:) .* sol.DT.Points);
    triplot(TR, '-k');
    
    for SolInd = [0 -1 1]
        SolPoint = sol.ApproxPoint(sol.ApproxPoint(:,end) == SolInd,1:end-1);
    
        scatter(ax, ...
            real(SolPoint), ...
            imag(SolPoint), ...
            25, ColorSet(SolInd+2,:), 'filled');
    end
    
    hold(ax, 'off')

    set(ax, FontSize=Global_FontSize);
    xlabel(ax, '$x$', FontSize=Axes_FontSize);
    ylabel(ax, '$y$', FontSize=Axes_FontSize);
    title(ax, Title_Set(i), FontSize=Title_FontSize);
    
    xlim(x);
    ylim(y);
end

%%

clear('sol')

%%

f = @(z) z + exp(4 * pi * z);

x = [-1 1];
y = x;

Tol = 1e-3;

%%

PointNumLim = [5e1 1e2 2e2 floor(1 / prod(Tol) / Tol(1))];

%%

fig = figure(Position=[0 0 1000 1000]);

t = 0;

for i = 1:4
    tic

    %%

    if exist('sol', 'var')
        sol.PointNumLim = PointNumLim(i);
        sol.fitTriang;
    else
        sol = GAMES(f, [x; y]', Tol, PointNumLim(i), Display=0);
    end

    %%

    t = t + toc;

    fprintf('Step <strong>#%i</strong> done!\nOverall time: <strong>%f sec</strong>\nPoints number: <strong>%i</strong>\n', i, t, sol.PointNum);

    %%

    ax = subplot(2,2,i);

    hold(ax, 'on')
    
    TR = triangulation(sol.DT(:,:), sol.GridNorm(1,:) + sol.GridNorm(2,:) .* sol.DT.Points);
    triplot(TR, '-k');
    
    for SolInd = [0 -1 1]
        SolPoint = sol.ApproxPoint(sol.ApproxPoint(:,end) == SolInd,1:end-1);
    
        scatter(ax, ...
            real(SolPoint), ...
            imag(SolPoint), ...
            25, ColorSet(SolInd+2,:), 'filled');
    end
    
    hold(ax, 'off')

    set(ax, FontSize=Global_FontSize);
    xlabel(ax, '$x$', FontSize=Axes_FontSize);
    ylabel(ax, '$y$', FontSize=Axes_FontSize);
    title(ax, Title_Set(i), FontSize=Title_FontSize);
    
    xlim(x);
    ylim(y);
end