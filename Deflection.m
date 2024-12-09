function [varargout] = Deflection(varargin)
%Deflection computes the deflection of a statically indeterminate beam.
%using finite element method.


N        = nargin;
Name     = varargin{1};
EI       = varargin{2};
Support  = varargin{3}; %location of support and type (1= clamped, 2= fixed);
Nodes    = varargin{4};

% Interpolating functions for deflection
N_1      = @(s) 1-3*s.^2 + 2*s.^3;
N_2      = @(s) s-2*s.^2 + s.^3;
N_3      = @(s) 3*s.^2 - 2*s.^3;
N_4      = @(s) -s.^2 + s.^3;
% Interpolating functions for slope
N_5      = @(s) -6*s + 6*s.^2;
N_6      = @(s) 1-4*s + 3*s.^2;
N_7      = @(s) 6*s - 6*s.^2;
N_8      = @(s) -2*s + 3*s.^2;

for n = 5:N
    LoadCell = varargin{n};
    type     = LoadCell{1};
    load     = LoadCell{2};
    loc      = LoadCell{3};
    if strcmp(type,'DF')
        if (numel(load) == 1) 
            load = repmat(load,size(loc)); 
        end
        p           = polyfit(loc, load, numel(loc) - 1); 
        varargin{n} = [varargin{n}, {p}];
    end
end
Kstiffness = zeros(2*numel(Nodes)); 
Fload      = zeros(2*numel(Nodes),1);
for i = 1:numel(Nodes)-1
    indx = (2*i - 1):(2*i + 2); 
    L    = Nodes(i+1) - Nodes(i);
    fci  = [0;0;0;0];
    K    = EI/L^3*[ 12    6*L    -12    6*L;
                    6*L   4*L^2  -6*L   2*L^2;
                   -12   -6*L    12    -6*L;
                    6*L   2*L^2  -6*L   4*L^2];
    for n = 5:N
        LoadCell = varargin{n};
        type     = LoadCell{1};
        load     = LoadCell{2};
        loc      = LoadCell{3};
        locgrid  = linspace(loc(1), loc(end),1000);
        if any((locgrid > Nodes(i)) & (locgrid <= Nodes(i+1)))
            if strcmp(type,'MF')
                s   = (loc - Nodes(i))/L;
                fci = fci+ [load*N_5(s)/L
                            load*N_6(s)
                            load*N_7(s)/L
                            load*N_8(s)];
            elseif strcmp(type,'CF')
                s   = (loc - Nodes(i))/L;
                fci = fci+ [load*N_1(s)
                           L*load*N_2(s)
                           load*N_3(s)
                           L*load*N_4(s)];
            elseif strcmp(type,'DF')
                start = max(Nodes(i),loc(1));
                stop  = min(Nodes(i+1),loc(end));
                if(start ~= stop)
                    s1    = (start - Nodes(i))/L;
                    s2    = (stop - Nodes(i))/L;
                    P     = LoadCell{4};
                    funP  = @(s)polyval(P,Nodes(i)+s*L);
                    fci   = fci + [L*integrator(funP,N_1, s1, s2)
                                   L^2*integrator(funP,N_2, s1, s2)
                                   L*integrator(funP,N_3, s1, s2)
                                   L^2*integrator(funP,N_4, s1, s2)];
                end
            end
        end
    end
    Fload(indx) = Fload(indx) + fci;
    Kstiffness(indx,indx) = Kstiffness(indx,indx) + K;
end
nodesolved = 1:2*numel(Nodes);
for i = numel(Nodes):-1:1
    node = Nodes(i);
    indx = find(node == Support(:,1));
    if(~isempty(indx))
        typ = Support(indx,2);
        if(typ == 1)
            deleteindex = 2*i-1;
        else
            deleteindex = [2*i-1,2*i];
        end
        Kstiffness(deleteindex,:) = [];
        Kstiffness(:,deleteindex) = [];
        Fload(deleteindex) = []; 
        nodesolved(deleteindex) = []; 
    end
end
vtheta = Kstiffness\Fload;
Sol             = zeros(2*numel(Nodes),1);
Sol(nodesolved) = vtheta;
%%
X = linspace(Nodes(1),Nodes(end),1000);
Y = zeros(size(X)); T = zeros(size(X)); 
v = Sol(1:2:end); t = Sol(2:2:end);
for i = 1:numel(Nodes)-1
    L  = Nodes(i+1) - Nodes(i);
    x  = X(Nodes(i) <= X & X <= Nodes(i+1)); 
    s  = (x - Nodes(i))/L;
    y  = N_1(s)*v(i) + L*N_2(s)*t(i) + N_3(s)*v(i+1) + L*N_4(s)*t(i+1);
    tt = N_5(s)*v(i)/L + N_6(s)*t(i) + N_7(s)*v(i+1)/L + N_8(s)*t(i+1);
    Y(Nodes(i) <= X & X <= Nodes(i+1)) = y;
    T(Nodes(i) <= X & X <= Nodes(i+1)) = tt; 
end
varargout = {v,t};
f = figure(1); f.Color = [1,1,1]; f.Position = [241 431 560 196];
plot(X,[Y;T],'linewidth',2); grid on;
legend('\delta(x)','\theta(x)', 'interpreter','latex', 'FontSize',12, ....
    'FontWeight','bold', 'location', 'best');
title(['$$\delta_{tip} = ', num2str(Y(end)),', \theta_{tip} = ',num2str(T(end)),...
          '$$'], 'interpreter','latex', 'FontSize',12, 'FontWeight','bold')
fig = getframe(gcf); imwrite(fig.cdata,[Name,'.png']);


function I = integrator(funP, funS, s1, s2)
[gls, glw] = Gaulegwt(s1, s2, 5);
I = 0;
for i = 1:5
    I = I + glw(i)*funP(gls(i))*funS(gls(i));
end