% This script is meant to explore other options for models of the Allee
% effect (specifically alternate forms of the weak model/other ones
% presented in literature)
close all; clear all; clc
%dN/dt = g(N-A)
%dN/dt = gN(N/A-1)
%dN/dt = gN-N(A/(1+BN))
%dN/dt = gN(1-(a+tau)/(N+tau))
%dN/dt = gN(1-(A/N)^p)

% Simulate each trajectory
%% Part 1: dN/dt = g(N-A)
tint = 1;
tsamp = 0:tint:200;
N1_init = [3; 8; 16];
g= 0.0238;
A =0;

f = @(t,N1) g*(N1-A);  % dN1/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1);

ttot = [];
N1tot = [];
for i = 1:length(N1_init)
[t,N1]=ode45(f, tsamp,N1_init(i), options);
ttot = horzcat(ttot,t);
N1tot = horzcat(N1tot,N1);

end

figure;
subplot(5,3,1)
for i = 1:length(N1_init)
plot(ttot(:,i), N1tot(:,i), 'c-', 'LineWidth',3)
text(ttot(55,i), N1tot(55, i), ['N_{0}=', num2str(N1_init(i))])
hold on
%text(ttot(end-20,i), N1tot(end-20,i),['Ni=', num2str(N1_init(i))]','HorizontalAlignment','left','VerticalAlignment','bottom','color','k' )
end
xlabel('time (hours)', 'FontSize', 16)
ylabel('Cell Number (N(t))', 'FontSize', 16)
ylim([ 0 300])
%title('dN/dt=gN solution')
title('Growth trajectory', 'FontSize',14)
subplot(5,3,3)
i=1
text(ttot(55,i)+5, log(N1tot(55, i)./(N1_init(i))), ['N_{0}=', num2str(N1_init(i)),', ', num2str(N1_init(i+1)),', &', num2str(N1_init(i+2))])
hold on
for i = 1:length(N1_init)
plot(ttot(:,i), log(N1tot(:,i)./(N1_init(i))), 'c-', 'LineWidth',3)
%text(ttot(55,i), log(N1tot(55, i)./(N1_init(i))), ['N_{0}=', num2str(N1_init(i)),',', num2str(N1_init(i+1),'&', num2str(N1_init(i+2))])
hold on
end
ylim([ 0 5])
xlabel('time (hours)', 'FontSize', 16)
ylabel('log(N(t)/N_{0})', 'FontSize', 16)
%title('Normalized growth rate dN/dt=g(N)')
title('Normalized growth rate', 'FontSize',14)
% Make a plot of instantaneous growth rate in time
for i = 1:length(N1_init)
    for j = 1:length(tsamp)-1
    percapg(j,i) = (N1tot(j+1,i)-N1tot(j,i))./(N1tot(j+1,i)*tint);
    end
end

subplot(5,3,2)
for i = 2%1:length(N1_init)
    %subplot(1,length(N1_init),i)
   plot(ttot(1:end-1,i), percapg(:,i), 'c-', 'LineWidth',3)
   text(ttot(55,i), percapg(75,i)+.001, ['N_{0}=',num2str(N1_init(i) )])
   hold off
   xlim([ 0 tsamp(end)])
   ylim([ 0 0.05])
   xlabel ('time (hours)', 'FontSize', 16)
   ylabel('^({^{\partialN}/_{\partialt}})/_{N}', 'FontSize',14)
   title(['Per capita growth rate'], 'FontSize', 14)
   %title(['Per capita growth rate for dN/dt =gN, N_{0}=',num2str(N1_init(i)),', A=', num2str(A)])
end


tint = 1;
tsamp = 0:tint:200;
N1_init = [3; 8; 16];
g= 0.0238;
A =5;

f = @(t,N1) g*(N1-A);  % dN1/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1);

ttot = [];
N1tot = [];
for i = 1:length(N1_init)
[t,N1]=ode45(f, tsamp,N1_init(i), options);
ttot = horzcat(ttot,t);
N1tot = horzcat(N1tot,N1);

end


subplot(5,3,4)
for i = 1:length(N1_init)
plot(ttot(:,i), N1tot(:,i), 'm-', 'LineWidth',3)
text(ttot(55,i), N1tot(55, i), ['N_{0}=', num2str(N1_init(i))])
hold on
%text(ttot(end-20,i), N1tot(end-20,i),['Ni=', num2str(N1_init(i))]','HorizontalAlignment','left','VerticalAlignment','bottom','color','k' )
end
xlabel('time (hours)', 'FontSize', 16)
ylabel('Cell Number (N(t))', 'FontSize', 16)
ylim([ 0 300])
%title('dN/dt=g(N solution')
%title('Strong Allee solution, A = 5')
subplot(5,3,6)
i=1
text(ttot(55,i)+5, log(N1tot(55, i)./(N1_init(i))), ['N_{0}=', num2str(N1_init(i)),', ', num2str(N1_init(i+1)),', &', num2str(N1_init(i+2))])
hold on
for i = 1:length(N1_init)
plot(ttot(:,i), log(N1tot(:,i)./(N1_init(i))), 'm-', 'LineWidth',3)
%text(ttot(55,i), log(N1tot(55, i)./(N1_init(i))), ['N_{0}=', num2str(N1_init(i))])
hold on
end
text(ttot(55,2), log(N1tot(55, 2)./(N1_init(2))), ['N_{0}=', num2str(N1_init(2))])
text(ttot(55,3), log(N1tot(55, 3)./(N1_init(2)))+0.3, ['N_{0}=', num2str(N1_init(3))])
ylim([ 0 5])
xlabel('time (hours)', 'FontSize', 16)
ylabel('log(N(t)/N_{0})', 'FontSize', 16)
%title('Normalized growth rate dN/dt=g(N)')
%title('Normalized growth rate dN/dt=g(N-A)')
% Make a plot of instantaneous growth rate in time
for i = 1:length(N1_init)
    for j = 1:length(tsamp)-1
    percapg(j,i) = (N1tot(j+1,i)-N1tot(j,i))./(N1tot(j+1,i)*tint);
    end
end

subplot(5,3,5)
for i = 2%1:length(N1_init)
    %subplot(1,length(N1_init),i)
   plot(ttot(1:end-1,i), percapg(:,i), 'm-', 'LineWidth',3)
   text(ttot(55,i), percapg(75,i)+.002, ['N_{0}=',num2str(N1_init(i) )])
   hold off
   xlim([ 0 tsamp(end)])
   ylim([ 0 0.05])
   xlabel ('time (hours)', 'FontSize', 16)
   ylabel('^{(^{\partialN}/_{\partialt})}/_{N}', 'FontSize',14)
   %title(['Per capita growth rate for dN/dt =g(N-A),A = 5'])
   %title(['Per capita growth rate for dN/dt =gN, N_{0}=',num2str(N1_init(i)),', A=', num2str(A)])
end

tint = 1;
tsamp = 0:tint:200;
N4_init = [3, 8, 16];
g= 0.0238;
a = -2;
tau = 6;


f = @(t,N4) g*N4*(1-((a+tau)./(N4+tau)));  % dN4/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1);

ttot = [];
N4tot = [];
for i = 1:length(N4_init)
[t,N4]=ode45(f, tsamp, N4_init(i), options);
ttot = horzcat(ttot,t);
N4tot = horzcat(N4tot, N4);

end

subplot(5,3,7)
for i = 1:length(N4_init)
plot(ttot(:,i), N4tot(:,i), 'g-', 'LineWidth',3)
hold on
text(ttot(55,i), N4tot(55, i), ['N_{0}=', num2str(N4_init(i))])
%text(ttot(end-20,i), N1tot(end-20,i),['Ni=', num2str(N1_init(i))]','HorizontalAlignment','left','VerticalAlignment','bottom','color','k' )
end
xlabel('time (hours)', 'FontSize', 16)
ylabel('Cell Number (N(t))', 'FontSize', 16)
ylim([ 0 300])
%title(['Weak Allee solution, A=', num2str(a),', \tau =', num2str(tau)])

subplot(5,3,9)
for i = 1:length(N4_init)
plot(ttot(:,i), log(N4tot(:,i)./(N4_init(i))), 'g-', 'LineWidth',3)
% text(ttot(55,i), log(N4tot(55, i)./(N4_init(i))), ['N_{0}=', num2str(N1_init(i))])
hold on
end
text(ttot(55,1), log(N4tot(55, 1)./(N4_init(1))), ['N_{0}=', num2str(N1_init(1))])
text(ttot(65,2), log(N4tot(65, 2)./(N4_init(2))), ['N_{0}=', num2str(N1_init(2))])
text(ttot(75,3), log(N4tot(75, 3)./(N4_init(3)))+.3, ['N_{0}=', num2str(N1_init(3))])
xlabel('time (hours)', 'FontSize', 16)
ylabel('log(N(t)/N_{0})', 'FontSize', 16)
%title(['Normalized growth rate dN/dt =gN(1-((A+\tau)/(N+\tau))), A=', num2str(a),', \tau =', num2str(tau)])
% Make a plot of instantaneous growth rate in time
for i = 1:length(N4_init)
    for j = 1:length(tsamp)-1
    percapg(j,i) = (N4tot(j+1,i)-N4tot(j,i))./(N4tot(j+1,i)*tint);
    end
end

subplot(5,3,8)
for i = 2%1:length(N1_init)
    %subplot(1,length(N1_init),i)
   plot(ttot(1:end-1,i), percapg(:,i), 'g-', 'LineWidth',3)
   hold on
   text(ttot(55,i), percapg(75,i)+.002, ['N_{0}=',num2str(N4_init(i) )])
   hold off
   xlim([ 0 tsamp(end)])
   ylim([ 0 0.05])
   xlabel ('time (hours)', 'FontSize', 16)
   ylabel('^{(^{\partialN}/_{\partialt})}/_{N}', 'FontSize',14)
   %title(['Per capita growth rate for dN/dt = gN(1-((A+\tau)/(N+\tau))), N_{0}=',num2str(N4_init(i)),', a=', num2str(a),', \tau =', num2str(tau)])
end


tint = 1;
tsamp = 0:tint:200;
N4_init = [3, 8, 16];
g= 0.0238;
a = -2;
tau = 6;
c = 4;

f = @(t,N4) g*N4*(1-((a+tau-c)./(N4+tau+c)));  % dN4/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1);

ttot = [];
N4tot = [];
for i = 1:length(N4_init)
[t,N4]=ode45(f, tsamp, N4_init(i), options);
ttot = horzcat(ttot,t);
N4tot = horzcat(N4tot, N4);

end

subplot(5,3,10)
for i = 1:length(N4_init)
plot(ttot(:,i), N4tot(:,i), 'k-', 'LineWidth',3)
hold on
text(ttot(55,i), N4tot(55, i), ['N_{0}=', num2str(N4_init(i))])
%text(ttot(end-20,i), N1tot(end-20,i),['Ni=', num2str(N1_init(i))]','HorizontalAlignment','left','VerticalAlignment','bottom','color','k' )
end
xlabel('time (hours)', 'FontSize', 16)
ylabel('Cell Number (N(t))', 'FontSize', 16)
ylim([ 0 300])
%title(['Weak Allee solution, A=', num2str(a),', \tau =', num2str(tau)])

subplot(5,3,12)
for i = 1:length(N4_init)
plot(ttot(:,i), log(N4tot(:,i)./(N4_init(i))), 'k-', 'LineWidth',3)
% text(ttot(55,i), log(N4tot(55, i)./(N4_init(i))), ['N_{0}=', num2str(N1_init(i))])
hold on
end
text(ttot(55,1), log(N4tot(55, 1)./(N4_init(1))), ['N_{0}=', num2str(N1_init(1))])
text(ttot(65,2), log(N4tot(65, 2)./(N4_init(2))), ['N_{0}=', num2str(N1_init(2))])
text(ttot(75,3), log(N4tot(75, 3)./(N4_init(3)))+.3, ['N_{0}=', num2str(N1_init(3))])
xlabel('time (hours)', 'FontSize', 16)
ylabel('log(N(t)/N_{0})', 'FontSize', 16)
%title(['Normalized growth rate dN/dt =gN(1-((A+\tau)/(N+\tau))), A=', num2str(a),', \tau =', num2str(tau)])
% Make a plot of instantaneous growth rate in time
for i = 1:length(N4_init)
    for j = 1:length(tsamp)-1
    percapg(j,i) = (N4tot(j+1,i)-N4tot(j,i))./(N4tot(j+1,i)*tint);
    end
end

subplot(5,3,11)
for i = 2%1:length(N1_init)
    %subplot(1,length(N1_init),i)
   plot(ttot(1:end-1,i), percapg(:,i), 'k-', 'LineWidth',3)
   hold on
   text(ttot(55,i), percapg(75,i)+.002, ['N_{0}=',num2str(N4_init(i) )])
   hold off
   xlim([ 0 tsamp(end)])
   ylim([ 0 0.05])
   xlabel ('time (hours)', 'FontSize', 16)
   ylabel('^{(^{\partialN}/_{\partialt})}/_{N}', 'FontSize',14)
   %title(['Per capita growth rate for dN/dt = gN(1-((A+\tau)/(N+\tau))), N_{0}=',num2str(N4_init(i)),', a=', num2str(a),', \tau =', num2str(tau)])
end

tint = 1;
tsamp = 0:tint:200;
N4_init = [3, 8, 16];
g= 0.0238;
a = -2;
tau = 6;
c = 4;

f = @(t,N4) g*N4*(1-((a+tau+c)./(N4+tau-c)));  % dN4/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1);

ttot = [];
N4tot = [];
for i = 1:length(N4_init)
[t,N4]=ode45(f, tsamp, N4_init(i), options);
ttot = horzcat(ttot,t);
N4tot = horzcat(N4tot, N4);

end

subplot(5,3,13)
for i = 1:length(N4_init)
plot(ttot(:,i), N4tot(:,i), 'r-', 'LineWidth',3)
hold on
text(ttot(55,i), N4tot(55, i), ['N_{0}=', num2str(N4_init(i))])
%text(ttot(end-20,i), N1tot(end-20,i),['Ni=', num2str(N1_init(i))]','HorizontalAlignment','left','VerticalAlignment','bottom','color','k' )
end
xlabel('time (hours)', 'FontSize', 16)
ylabel('Cell Number (N(t))', 'FontSize', 16)
ylim([ 0 300])
%title(['Weak Allee solution, A=', num2str(a),', \tau =', num2str(tau)])

subplot(5,3,15)
for i = 1:length(N4_init)
plot(ttot(:,i), log(N4tot(:,i)./(N4_init(i))), 'r-', 'LineWidth',3)
% text(ttot(55,i), log(N4tot(55, i)./(N4_init(i))), ['N_{0}=', num2str(N1_init(i))])
hold on
end
text(ttot(55,1), log(N4tot(55, 1)./(N4_init(1))), ['N_{0}=', num2str(N1_init(1))])
text(ttot(65,2), log(N4tot(65, 2)./(N4_init(2))), ['N_{0}=', num2str(N1_init(2))])
text(ttot(75,3), log(N4tot(75, 3)./(N4_init(3)))+.3, ['N_{0}=', num2str(N1_init(3))])
xlabel('time (hours)', 'FontSize', 16)
ylabel('log(N(t)/N_{0})', 'FontSize', 16)
%title(['Normalized growth rate dN/dt =gN(1-((A+\tau)/(N+\tau))), A=', num2str(a),', \tau =', num2str(tau)])
% Make a plot of instantaneous growth rate in time
for i = 1:length(N4_init)
    for j = 1:length(tsamp)-1
    percapg(j,i) = (N4tot(j+1,i)-N4tot(j,i))./(N4tot(j+1,i)*tint);
    end
end

subplot(5,3,14)
for i = 2%1:length(N1_init)
    %subplot(1,length(N1_init),i)
   plot(ttot(1:end-1,i), percapg(:,i), 'r-', 'LineWidth',3)
   hold on
   text(ttot(55,i), percapg(75,i)+.002, ['N_{0}=',num2str(N4_init(i) )])
   hold off
   xlim([ 0 tsamp(end)])
   ylim([ 0 0.05])
   xlabel ('time (hours)', 'FontSize', 16)
   ylabel('^{(^{\partialN}/_{\partialt})}/_{N}', 'FontSize',14)
   %title(['Per capita growth rate for dN/dt = gN(1-((A+\tau)/(N+\tau))), N_{0}=',num2str(N4_init(i)),', a=', num2str(a),', \tau =', num2str(tau)])
end


