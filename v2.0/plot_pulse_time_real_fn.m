function plot_pulse_time_real_fn(net,f_ind)

if net(1).param.re_var_max>0
    figure
    subplot(2,length(net(f_ind).s),1)
    title(strcat('freq=',num2str(net(1).f_var(f_ind))))
    hold on;grid on;
    plot(net(f_ind).t_var,real(net(f_ind).curr.Ve),'linewidth',2);
    text(net(f_ind).t_var(end/2),net(f_ind).curr.Vdc_e_r,strcat('Vdc=',num2str(net(f_ind).curr.Vdc_e_r),'V, dVpp=',num2str(net(f_ind).curr.dVpp_e_r*1e6),'uV, PAT=',num2str(net(f_ind).curr.PAT_e_r*1e3),'ms'),'FontSize',13);
    xlabel('time');ylabel('Re(I_{Ve}) (V)')
    ax = gca; ax.FontSize=13;
    
    
    for i=1:length(net(f_ind).s)
        subplot(2,length(net(f_ind).s),length(net(f_ind).s)+i)
        hold on;grid on;
        plot(net(f_ind).t_var,real(net(f_ind).s(i).Ve)*1e3,'linewidth',2);
        text(net(f_ind).t_var(end/2),net(f_ind).s(i).Vdc_e_r*1e3,strcat('Vdc=',num2str(net(f_ind).s(i).Vdc_e_r*1e3),'mV, dVpp=',num2str(net(f_ind).s(i).dVpp_e_r*1e6),'uV, PAT=',num2str(net(f_ind).s(i).PAT_e_r*1e3),'ms'),'FontSize',13);
        xlabel('time');ylabel('Re(V_{Ve}) (mV)')
        ax = gca; ax.FontSize=13;
    end
    
    
    if net(1).param.write_figures==1 && net(1).param.debug==0
        file_str=strcat(net(1).param.path);
        print(strcat(file_str,'pulse_time_real.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
        %savefig(strcat(file_str,'pulse_time_real.fig')); %<-Save as PNG with 300 DP
    end
end