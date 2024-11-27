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

f = @(z, p, m, n, T) D(z .* exp(1i * pi / 4), p, m, n, T);

x = [1e-1 2];
y = [-5 0];
p = [0 2];

Tol = [1e-3 1e-2];

%%

T = [
    1 1 1;
    1 10 1
    ];
PointNumLim = [1e6 1.5e6];

%%

fig = figure(Position=[0 0 1000 1000]);

for i = 1:2
    tic

    %%

    sol = GAMES(@(v) f(v(:,1), v(:,2), m, n, T(i,:)), [x; y; p]', Tol, PointNumLim(i));

    %%

    fprintf('Step <strong>#%i</strong> done!\nOverall time: <strong>%f sec</strong>\nPoints number: <strong>%i</strong>\n', i, toc, sol.PointNum);

    %%

    ax = subplot(2,2,i+mod(i-1,2));

    hold(ax, 'on')
    
    for SolInd = [0 -1 1]
        SolPoint = sol.ApproxPoint(sol.ApproxPoint(:,end) == SolInd,1:end-1);
        SolPoint(:,1) = SolPoint(:,1) .* exp(1i * pi / 4) / sqrt(2);
    
        scatter3(ax, ...
            real(SolPoint(:,1)), ...
            imag(SolPoint(:,1)), ...
            SolPoint(:,2), ...
            5, ColorSet(SolInd+2,:), 'filled');
    end
    
    hold(ax, 'off')

    set(ax, FontSize=Global_FontSize);
    xlabel(ax, '$\Omega/\omega_{pe}$', FontSize=Axes_FontSize);
    ylabel(ax, '$\gamma/\omega_{pe}$', FontSize=Axes_FontSize);
    zlabel(ax, '$k/k_{Dc}$', FontSize=Axes_FontSize);
    title(ax, Title_Set(i+mod(i-1,2)), FontSize=Title_FontSize);

    view(ax,0,0);
    
    xlim([0 3.5 / sqrt(2)]);
    % ylim([-3 / sqrt(2) 0]);
    zlim([0 2]);

    %%

    ax = subplot(2,2,i+mod(i-1,2)+1);

    hold(ax, 'on')
    
    for SolInd = [-1 0 1]
        SolPoint = sol.ApproxPoint(sol.ApproxPoint(:,end) == SolInd,1:end-1);
        SolPoint(:,1) = SolPoint(:,1) .* exp(1i * pi / 4) / sqrt(2);
    
        scatter3(ax, ...
            real(SolPoint(:,1)), ...
            imag(SolPoint(:,1)), ...
            SolPoint(:,2), ...
            5, ColorSet(SolInd+2,:), 'filled');
    end
    
    hold(ax, 'off')

    set(ax, FontSize=Global_FontSize);
    set(ax, 'YDir', 'reverse');
    xlabel(ax, '$\Omega/\omega_{pe}$', FontSize=Axes_FontSize);
    ylabel(ax, '$\gamma/\omega_{pe}$', FontSize=Axes_FontSize);
    zlabel(ax, '$k/k_{Dc}$', FontSize=Axes_FontSize);
    title(ax, Title_Set(i+mod(i-1,2)+1), FontSize=Title_FontSize);

    view(ax,90,0);
    
    % xlim([0 3.5 / sqrt(2)]);
    ylim([-3 / sqrt(2) 0]);
    zlim([0 2]);
end

%%

function value = D(omega, k, m, n, T)
    m = reshape(m,1,1,[]);
    n = reshape(n,1,1,[]);
    T = reshape(T,1,1,[]);
    omegaZero = (omega == 0);
    kZero = (k == 0);

    value = zeros(size(omega));
    value(~omegaZero & kZero) = 1 - sum(n ./ m,3) ./ omega(~omegaZero & kZero).^2;
    value(~kZero) = 1 - 1 ./ (2 * k(~kZero).^2) .* sum(n ./ T .* pdfZDiv(omega(~kZero) ./ k(~kZero) .* sqrt(m ./ (2 * T))),3);
end