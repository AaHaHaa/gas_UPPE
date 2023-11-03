function pcolor_plotbox(ax)
%PCOLOR_PLOTBOX It adds box lines to the pcolor figure
%   ax: the axis handle of the pcolor figure

Xlim = ax.XLim;
Ylim = ax.YLim;
hold on;
hx1 = plot(linspace(Xlim(1),Xlim(2),10),Ylim(1)*ones(1,10)); set(hx1,'linewidth',1,'color','k');
hx2 = plot(linspace(Xlim(1),Xlim(2),10),Ylim(2)*ones(1,10)); set(hx2,'linewidth',1,'color','k');
hy1 = plot(Xlim(1)*ones(1,10),linspace(Ylim(1),Ylim(2),10)); set(hy1,'linewidth',1,'color','k');
hy2 = plot(Xlim(2)*ones(1,10),linspace(Ylim(1),Ylim(2),10)); set(hy2,'linewidth',1,'color','k');
hold off;

end

