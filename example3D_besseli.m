clear
set(0, defaulttextinterpreter='latex')
addpath('specialFunc')

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

m = [1 1 1836];
n = [1 1 2];

x = [1 6];
y = [-0.5 0.5];
p = [1e-2 10];

Tol = [1e-3 1e-2];
PointNumLim = 2.5e4;

%%

p2c = [2 4];

%%

fig = figure(Position=[0 0 1000 1000]);

for i = 1:2
    Res = [];
    TotalPoint = 0;

    tic

    %%

    for n = min(x):max(x)-1
        sol = GAMES(@(v) D(v(:,1), v(:,2), p2c(i), 10), [n n + 1; y; p]', Tol, PointNumLim, StopCond=0);
        Res = ...
            [
            Res;
            sol.ApproxPoint
            ];
        TotalPoint = TotalPoint + sol.PointNum;
    end

    %%

    fprintf('Step <strong>#%i</strong> done!\nOverall time: <strong>%f sec</strong>\nPoints number: <strong>%i</strong>\n', i, toc, TotalPoint);

    %%

    ax = subplot(2,1,i);

    hold(ax, 'on')
    
    for SolInd = [0 -1 1]
        SolPoint = Res(Res(:,end) == SolInd,1:end-1);
    
        scatter3(ax, ...
            real(SolPoint(:,1)), ...
            imag(SolPoint(:,1)), ...
            SolPoint(:,2), ...
            5, ColorSet(SolInd+2,:), 'filled');
    end
    
    hold(ax, 'off')

    set(ax, FontSize=Global_FontSize);
    xlabel(ax, '$\Omega/\omega_c$', FontSize=Axes_FontSize);
    ylabel(ax, '$\gamma/\omega_c$', FontSize=Axes_FontSize);
    zlabel(ax, '$kv/\omega_c$', FontSize=Axes_FontSize);
    title(ax, Title_Set(i), FontSize=Title_FontSize);

    view(ax,0,0);
    
    xlim(x);
    ylim(y);
    zlim(p);
end

%%

function value = D(omega, k, p2c, n_max)
    value = 1 - besseli(0, k.^2) .* exp(-k.^2);

    for n=1:n_max
        value = value - (besseli(n, k.^2) ./ (omega - n) + besseli(-n, k.^2) ./ (omega + n)) .* omega .* exp(-k.^2);
    end

    value = 1 + p2c^2 ./ k.^2 .* value;
end