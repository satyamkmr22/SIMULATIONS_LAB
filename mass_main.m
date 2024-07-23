clc;
clear all;

% Define the data
xdata = [0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.0];
ydata = [0, 0.134, 0.23, 0.304, 0.365, 0.418, 0.579, 0.665, 0.729, 0.779, 0.825, 0.87, 0.915, 0.958, 0.979, 1];
x_points = [0.97, 0.45];
y_points = [0.97, 0.45];

% Define the function
fun = @(param, xdata) (param(1) .* xdata) ./ (1 + (param(2) .* xdata) + (param(3) .* (xdata.^2)));

% Initial guess for parameters
initial_guess = [1, 1, 1];

% Perform curve fitting
param = lsqcurvefit(fun, initial_guess, xdata, ydata);

% Extract parameters
a = param(1);
b = param(2);
c = param(3);
y = a .* xdata ./ (1 + b .* xdata + c .* (xdata.^2));

figure(1)
plot(xdata, ydata);
hold on;
plot(xdata, y);
xlabel('xdata');
ylabel('y')
legend('given', 'predicted');

xdata = linspace(0, 1, 200);
y_calc = (a .* xdata) ./ (1 + (b .* xdata) + (c .* (xdata.^2)));
y_normal = xdata;

figure(2)
plot(xdata, y_calc);
hold on;
plot(xdata, y_normal);
hold on;

% Plot the points
plot(x_points, y_points, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
hold on;

q = 0.8;
zf = 0.45;
y_pinch_point=a*0.7/(1+(b*0.7)+(c*0.7*0.7));
y_pinch_point
m_reflux=(0.97-y_pinch_point)/(0.97-0.7);
c_reflux=y_pinch_point-(m_reflux*0.7);
Reflux_hypo=(0.97/c_reflux)-1;
Reflux_actual=2.5*Reflux_hypo;

% Define the function to solve for the intersection of feed line and equilibrium curve
fun = @(x) a * x ./ (1 + (b * x) + c * (x.^2)) + (q / (1 - q)) * x - zf / (1 - q);
x0 = 0;
xdata_feed = linspace(0.37, 0.5, 50);
x_solution = fsolve(fun, x0);
y_feed = ((q / (q - 1)) * xdata_feed) - zf / (q - 1);
plot(xdata_feed, y_feed);
% Calculate y using the obtained x_solution
y_solution = a * x_solution / (1 + b * x_solution + c * x_solution.^2);
plot(x_solution, y_solution, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

% Calculate the equation of the line passing through (0.97, 0.97) and (x_solution, y_solution)
m = Reflux_actual/(Reflux_actual+1);
c_re = 0.97/(Reflux_actual+1);

% Draw the line
x_line = linspace(0.97, 0, 200);
y_line = m * x_line + c_re;
plot(x_line, y_line, 'k--');
hold on;


% Solve the system of equations
syms x y; % Define symbolic variables
eq1 = x + y  == 450;
eq2 = x * 0.97 + y * 0.02 + 35 == 225;
sol = solve([eq1, eq2], [x, y]);

% Extract solutions
D = sol.x; % flow rate of distillate
W = sol.y; % flow rate of bottom output
D=double(D);
W=double(W);

y_seventyfive_percent=m*0.70+c_re;
plot(0.70,y_seventyfive_percent,'o', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
y_seventyfive_percent
hold on;

x_line=linspace(0.97,0.7,100);
y_line = m * x_line + c_re;
plot(x_line,y_line,'k');
hold on;
y_line_side=linspace(0.7,y_pinch_point,80);
x_line_side=zeros(1,80);
x_line_side(1,:)=0.7;
plot(x_line_side, y_line_side, 'color', [0.5 0 0.5]);
hold on;
plot(0.70,0.7,'o', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
hold on;
plot(0.70,y_pinch_point,'o', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
hold on;
V=(Reflux_actual+1)*D;
S=50;
xs=0.70;
xd=0.97;
xw=0.02;
L=D*Reflux_actual;
L_dash=L-S;

% Define the equations
L_dash_over_V = L_dash / V;
numerator = D * xd + S * xs;
denominator = D * (Reflux_actual + 1);

% Define the line equation
m_section2=L_dash_over_V;
c_section2=numerator/denominator;
fun=@(x_feed) L_dash_over_V*x_feed+(numerator/denominator)+(4*x_feed)-(45/20);
x_feed=fsolve(fun,0.5);
y_feed_intersec=-(4*x_feed)+(45/20);
plot(x_feed,y_feed_intersec,'o', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

% Plot the line joining (0.70, 0.854) and (x_feed, y_feed_intersec)
plot([0.70, x_feed], [y_seventyfive_percent, y_feed_intersec], 'k');
plot([0.02, x_feed], [0.02, y_feed_intersec], 'k');
hold on;
plot(0.02,0.02,'o', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
hold on;

m_section3=(y_feed_intersec-0.02)/(x_feed-0.02);
c_section3=y_feed_intersec-(m_section3*x_feed);
x_for_tray=0.7;
x_values=zeros(1000,1);
count=0;
x_vary=0.70;
y_for_tray=y_seventyfive_percent;
y_vary=y_for_tray;
check=0;
count2=0;
while x_vary>0.02
eqb=@(x) y_vary-((a*x)/(1+(b*x)+(c*(x^2))));
x=fsolve(eqb,0.5);
x_prev=x_vary;
x_vary=x;
y_prev=y_vary;
if x_vary>x_feed
y_vary=(m_section2*x_vary)+c_section2;
else
y_vary=(m_section3*x_vary)+c_section3;
check=1;
end
if check==1
count2=count;
end
count=count+1;
plot([x_prev,x_vary],[y_prev,y_prev],'r');
hold on;
plot([x_vary,x_vary],[y_prev,y_vary],'b');
hold on;
legend('eqb curve','y=x','distillate and feed points','feed line','intersec feed & eqb curve','section 2 operating line','intersec side & section1','section 1 operating line','side line','side stream mole fraction','pinch point','section 2 intersec with feed','section 2 operating line','section 3 operating line','0.02,0.02');
xlabel('x');
ylabel('y');
end
x_final_store=x_vary;
x_final_pre=x_prev;
count
x_vary=0.7;
y_vary=y_seventyfive_percent;
count1=0;
while x_vary<0.97
    x_prev=x_vary;
    y_prev=y_vary;
y_vary=a*x_vary/(1+(b*x_vary)+(c*(x_vary^2)));
x_vary=(y_vary-c_re)/m;
count1=count1+1;
 % count1 is the side stream tray count
plot([x_prev,x_prev],[y_prev,y_vary],'r');
hold on;
plot([x_prev,x_vary],[y_vary,y_vary],'b');
hold on;
end
count=count-1;
count=count+(x_final_pre-0.02)/(x_final_pre-x_final_store)+count1;
count
count2=count1+count2;
count1
%feed tray is around 6 th
% side stream tray is 4th
%number of trays is 9.819821394239320
% L at section three boundary which is lying above the feed stage is 2.696620996396439e+02
% L at section two boundary which is lying above side stream stage is 2.196620996396439e+02
% V at both the sections is 4.601884154291176e+02 (section two and three)
L_below_feed=(0.8*500)+L_dash;
V_below_feed=V-(0.2*500);