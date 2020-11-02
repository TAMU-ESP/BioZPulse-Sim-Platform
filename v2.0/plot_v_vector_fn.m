function plot_v_vector_fn(net_arr)

name_base={'Object Loc. '};

% Plot Vs vector
i=1;
name_arr{i}=strcat(num2str([round(net_arr{i}(1).Ni)]'));
for j=1:4
    subplot(4,1,j)
    s_arr=cat(1,net_arr{i}.s);
    s_arr=permute(s_arr,[2,1]);
    Vs_arr{i}=[s_arr.Vs];
    x_axis{i}=(1:length(Vs_arr{i}));
    L=length(x_axis{i});
    h(i)=semilogy(x_axis{i}((j-1)*L/4+[1:L/4]),Vs_arr{i}((j-1)*L/4+[1:L/4]),'linewidth',2);
    grid on; hold on;
    Vs_max(i)=max(Vs_arr{i});
    Vs_min(i)=min(Vs_arr{i});
    hold on;
    if j==1
        title(strcat('Current Injection at side Y=0'));
    elseif j==2
        title(strcat('Current Injection at side X=',num2str(net_arr{i}(1,1).Nj)));
    elseif j==3
        title(strcat('Current Injection at side Y=',num2str(net_arr{i}(1,1).Ni)));
    else
        title(strcat('Current Injection at side X=0'));
    end
    for k=1:net_arr{i}(1,1).param.Nelec/4
        h(length(net_arr)+1)=semilogy(x_axis{i}((j-1)*L/4+1)-1+[net_arr{i}(1,1).param.Nelec-3 net_arr{i}(1,1).param.Nelec-3]*k,[min(Vs_min) max(Vs_max)],'--k');
    end
    legend(h(1:i+1),{name_arr{:},'Last voltage measurement\newline at a current injection location'},'fontsize',11,'Location','eastoutside')
    xlim([x_axis{i}((j-1)*L/4+1) x_axis{i}((j-1)*L/4+L/4)])
    if j==2
        ylabel('Electrode Voltage Difference Ve (mV)');
    end
    ax = gca; ax.FontSize=12;
    set(gca, 'FontName', 'Times New Roman');
end
xlabel('Voltage Measurement Index');


