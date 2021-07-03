%% Plotting
%Examples of cool plotting technqiues

%% Subplot

%plots
Fig=figure;
Ax1=subplot(2,1,1,'parent',Fig);
Ax2=subplot(2,1,2,'parsent',Fig);

hold(Ax1,'on'); hold(Ax2,'on'); 
grid(Ax1,'on'); grid(Ax2,'on');

xlabel(Ax1,'Time [s]'); ylabel(Ax1,'regular sine [m]');
xlabel(Ax2,'Time [s]'); ylabel(Ax2,'smoothed sine [m]');
title(Ax1,'first title'); title(Ax2,'second title');

t=linspace(0,10*pi,1000);
y1=sin(t);
y2=smooth(sin(t));

plot(Ax1,t,y1);  plot(Ax2,t,y2);

legend(Ax1,'blue line'); legend(Ax2,'blue line too');