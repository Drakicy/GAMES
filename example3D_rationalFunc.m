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

f = @(v) (v(:,1) + v(:,2) - (sqrt(3) + 1i * sqrt(2)) / 4) ./ (v(:,1) - v(:,2) + (sqrt(3) + 1i * sqrt(2)) / 4);

x = [-1 1];
y = x;
p = x;

Tol = [1e-3 1e-2];

%%

PointNumLim = [1e3 5e3 1e4 floor(1 / prod(Tol) / Tol(1))];

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
        sol = GAMES(f, [x; y; p]', Tol, PointNumLim(i), Display=0);
    end

    %%

    t = t + toc;

    fprintf('Step <strong>#%i</strong> done!\nOverall time: <strong>%f sec</strong>\nPoints number: <strong>%i</strong>\n', i, t, sol.PointNum);

    %%

    ax = subplot(2,2,i);

    hold(ax, 'on')
    
    for SolInd = [0 -1 1]
        SolPoint = sol.ApproxPoint(sol.ApproxPoint(:,end) == SolInd,1:end-1);
    
        scatter3(ax, ...
            real(SolPoint(:,1)), ...
            imag(SolPoint(:,1)), ...
            SolPoint(:,2), ...
            5, ColorSet(SolInd+2,:), 'filled');
    end

    hold(ax, 'off')

    set(ax, FontSize=Global_FontSize);
    xlabel(ax, '$x$', FontSize=Axes_FontSize);
    ylabel(ax, '$y$', FontSize=Axes_FontSize);
    zlabel(ax, '$p$', FontSize=Axes_FontSize);
    title(ax, Title_Set(i), FontSize=Title_FontSize);

    view(ax,45,30);
    
    xlim(x);
    ylim(y);
    zlim(p);
end

%%

clear('sol')

f = @(v) 4 * v(:,1).^2 + 4 * v(:,2).^2 - 1;

x = [-1 1];
y = x;
p = x;

Tol = [1e-3 1e-2];

%%

PointNumLim = [1e3 5e3 1e4 floor(1 / prod(Tol) / Tol(1))];

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
        sol = GAMES(f, [x; y; p]', Tol, PointNumLim(i), Display=0);
    end

    %%

    t = t + toc;

    fprintf('Step <strong>#%i</strong> done!\nOverall time: <strong>%f sec</strong>\nPoints number: <strong>%i</strong>\n', i, t, sol.PointNum);

    %%

    ax = subplot(2,2,i);

    hold(ax, 'on')
    
    for SolInd = [-1 0 1]
        SolPoint = sol.ApproxPoint(sol.ApproxPoint(:,end) == SolInd,1:end-1);
    
        scatter3(ax, ...
            real(SolPoint(:,1)), ...
            imag(SolPoint(:,1)), ...
            SolPoint(:,2), ...
            5, ColorSet(SolInd+2,:), 'filled');
    end

    hold(ax, 'off')

    set(ax, FontSize=Global_FontSize);
    xlabel(ax, '$x$', FontSize=Axes_FontSize);
    ylabel(ax, '$y$', FontSize=Axes_FontSize);
    zlabel(ax, '$p$', FontSize=Axes_FontSize);
    title(ax, Title_Set(i), FontSize=Title_FontSize);

    view(ax,45,30);
    
    xlim(x);
    ylim(y);
    zlim(p);
end