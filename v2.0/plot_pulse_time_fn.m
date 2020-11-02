function plot_pulse_time_fn(net,f_ind)

if net(1).param.re_var_max>0
    figure       
    subplot(2,2,1)
    for i=1:length(net(f_ind).s)       
        x=net(f_ind).t_var;
        y=(abs(net(f_ind).s(i).Ve)-abs(mean(net(f_ind).s(i).Ve)));
        xfit=[x(1):0.01:x(end)];
        p=spline(x,y);
        yfit=ppval(p,xfit);
    
        hold on;grid on;
        plot(xfit,yfit*1e6,'linewidth',2);
        %str=sprintf('{\it{V}%g}_{dc}=',i);
        %text(1.3,(3-i)*6,str,'FontSize',13, 'FontName', 'Times');
        str1=sprintf('%1.2f',net(f_ind).s(i).Vdc_e*1e3);
        str2=sprintf('%1.2f',net(f_ind).s(i).dVpp_e*1e6);
        str3=sprintf('%1.2f',net(f_ind).s(i).PAT_e*1e3);
        text(1.3,(3-i)*40,strcat('{\it{V}',num2str(i),'}_{dc}=',str1,'mV, {\it{\DeltaV}',num2str(i),'}_{pp}=',str2,'\muV, {\it{T}',num2str(i),'}_{d}=',str3,'ms'),'FontSize',13, 'FontName', 'Times');      
        %text(1.3,(3-i)*6,strcat('{\it{V}',num2str(i),'}_{dc}=',num2str(net(f_ind).s(i).Vdc_e*1e3),'mV, {\it{\DeltaV}',num2str(i),'}_{pp}=',num2str(net(f_ind).s(i).dVpp_e*1e6),'\muV, {\it{T}',num2str(i),'}_{d}=',num2str(net(f_ind).s(i).PAT_e*1e3),'ms'),'FontSize',13, 'FontName', 'Times');
        xlabel('Time (s)');ylabel('{\it{V}}_{Sim}(\it{t}) - {\it{V}}_{dc} (\muV)')
        legend('{\it{V1}}','{\it{V2}}','location','southeast')
        ax = gca; ax.FontSize=13;
        set(gca, 'FontName', 'Times New Roman');
    end
    set(gcf,'position',[0.1634    0.2078    1.2124    0.4640]*1e3)
    
    if net(1).param.write_figures==1 && net(1).param.debug==0
        file_str=strcat(net(1).param.path);
        print(strcat(file_str,'pulse_time.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
        %savefig(strcat(file_str,'pulse_time.fig')); %<-Save as PNG with 300 DP
    end
end