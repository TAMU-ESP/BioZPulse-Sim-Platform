%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%   
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation 
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling, 
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   Simulate DC Bio-Z(V_DC) and arterial pulse amplitude(deltaV_PP) versus frequency

clc
clear
close all
fclose('all');
warning('off','all')

write_figures=0;
write_mat=1;
write_netlist='freq';
% Dimensions
Nk=7;
% Current Source
i_num=1;
i_opp=0;
i_freq=100;
dia_arr=[0];
dia_ar=1;
s_num=1;
margin_x=15;
margin_y=6;
spacing_x=max(dia_arr)*dia_ar+8;
spacing_y=max(dia_arr)+[8];
%spacing_y_sns_arr=max(dia_arr)+[8 16 24];
spacing_y_sns_arr=max(dia_arr)+[10 20 30];
dia_x=dia_arr*dia_ar;
dia_y=dia_arr;

Nj=(max(dia_arr)*dia_ar+2*margin_x);
i_src_n1j_arr=round(Nj/2-dia_x/2);
% Configurations
z_en=1;
nstep=8;
mode=1;  %0=> .TRAN, 1=> .AC
ac_freq=[20 1.5e3 100e3]; %npts fstart fend
fi_req=[1,4];
debug=0;
sim_speed=3;
dist_step_count=10;
res=2;

for t1=1:length(spacing_y_sns_arr)
    clear net_arr2 net;
    spacing_y_sns=spacing_y_sns_arr(t1);
    Ni=(max(dia_arr)+(s_num*2-1)*spacing_y_sns+2*spacing_y+2*margin_y);
    spacing_arr_in=spacing_y;
    i_src_n1j=i_src_n1j_arr;
    % Sensing
    i_dia_x=dia_x;
    s_dia_x=dia_x;
    i_dia_y=dia_y;
    s_dia_y=dia_y;
    i_src_n1i=round((Ni-(dia_y+(s_num*2-1)*spacing_y_sns+2*spacing_y))/2);
    s_n1i=[];s_n1j=[];s_n2i=[];s_n2j=[];
    if s_num>0
        for i=1:s_num
            if s_num==1
                s_n1i(i)=i_src_n1i+spacing_y*(2*i-1);
                s_n2i(i)=s_n1i(i)+spacing_y_sns;
                s_n1j(i)=i_src_n1j;
                s_n2j(i)=s_n1j(i);
            else
                s_n1i(i)=s_n2i(i-1)+spacing_y*(2*i-1);
                s_n2i(i)=s_n1i(i)+spacing_y_sns;
                s_n1j(i)=i_src_n1j;
                s_n2j(i)=s_n1j(i);
            end
        end
    end
    i_src_n2i=s_n2i(end)+spacing_y;
    i_src_n2j=i_src_n1j;
    % Artery Model
    a=3;
    r_dc=30*a;
    ci_dc=100e-9/6*a;
    ri_dc=r_dc;   
    r_art_dc=4/7*r_dc;
    ri_art_dc=4/7*ri_dc;
    ci_art_dc=1/(4/7)*ci_dc;
    r_var_freq_arr=[1];
    r_var_time_arr=1./min(r_var_freq_arr,[],2)*1;
    S=size(r_var_freq_arr,2);
    net.param.S=S;
    r_var_delay_arr=50e-3*ones(1,S);
    r_var_n1j=[round((Nj)*(1/2))]-1;    % R var start location vertical
    r_var_n1i=1*ones(1,S)-1;    % R var start location horizontal
    r_var_n1k_arr=[3]*ones(1,S);     % R var start location depth
    r_var_max=r_art_dc*1e-2*ones(1,S);
    ri_var_max=ri_art_dc*1e-2*ones(1,S);
    ci_var_max=ci_art_dc*-1e-2*ones(1,S);
    r_var_steps=3*ones(1,S);
    r_var_v=0*ones(1,S);       % 1-> vertical, 0-> horizontal
    r_var_dim=3*ones(1,S);     % Dimension (1,2,3)
    r_var_nan=1*ones(1,S);     % 0-> No NaN, 1-> NaN
    r_var_dia=2*ones(1,S);

    if mode==0
        i_val=1e-3;
        i_ac=0;
        time_arr=0;
    else
        i_val=0;
        i_ac=1e-3;
        if debug==0
            time_arr=[0:r_var_time_arr/nstep:r_var_time_arr/nstep*(nstep-1)]';
        else
            time_arr=0;
        end
    end
    color_arr=[0    0.4470    0.7410;
        0.8500    0.3250    0.0980;
        0.9290    0.6940    0.1250;
        0.4940    0.1840    0.5560;
        0.4660    0.6740    0.1880;
        0.3010    0.7450    0.9330;
        0.6350    0.0780    0.1840;
        1.0000         0         0;
        0.2500    0.2500    0.2500;
        0    0.5000         0];
    
    tindex=0;
    for t0=1:length(time_arr)
        %% Generate and Simualte Spice
        
        r_var_delay=r_var_delay_arr(1);
        r_var_freq=r_var_freq_arr(1,:);
        r_var_time=r_var_time_arr(1);
        tstep=r_var_time*nstep;
        r_var_n1k=r_var_n1k_arr(1);
        time=time_arr(t0);
        
        net.Ni=round(Ni/res);
        net.Nj=round(Nj/res);
        net.Nk=round(Nk/res);
        net.param.i_src_n1i=round(i_src_n1i/res);
        net.param.i_src_n2i=round(i_src_n2i/res);
        net.param.i_src_n1j=round(i_src_n1j/res);
        net.param.i_src_n2j=round(i_src_n1j/res);
        net.param.i_num=i_num;
        net.param.i_opp=i_opp;
        net.param.r_dc=r_dc;
        net.param.ri_dc=ri_dc;
        net.param.ci_dc=ci_dc;
        net.param.r_art_dc=r_art_dc;
        net.param.ri_art_dc=ri_art_dc;
        net.param.ci_art_dc=ci_art_dc;
        net.param.r_var_max=r_var_max;
        net.param.ri_var_max=ri_var_max;
        net.param.ci_var_max=ci_var_max;
        net.param.r_var_steps=r_var_steps;
        net.param.r_var_n1i=round(r_var_n1i/res);
        net.param.r_var_n1j=round(r_var_n1j/res);
        net.param.r_var_n1k=round(r_var_n1k/res);
        net.param.r_var_v=r_var_v;
        net.param.r_var_dim=r_var_dim;
        net.param.r_var_num=net.Ni*ones(1,S);
        net.param.r_var_time=r_var_time;
        net.param.r_var_seg=net.Ni*ones(1,S);
        net.param.r_var_dia=round(r_var_dia/res);
        net.param.r_var_delay=r_var_delay;
        net.param.r_var_freq=r_var_freq;
        net.param.tstep=tstep;
        net.param.time=time;       
        net.curr(1).dia_x=round(i_dia_x/res);
        net.curr(1).dia_y=round(i_dia_y/res);
        net.curr(1).ival=i_val;
        net.curr(1).iac=i_ac;
        net.curr(1).freq=i_freq;        
        net.param.s_n1i=round(s_n1i/res);
        net.param.s_n2i=round(s_n2i/res);
        net.param.s_n1j=round(s_n1j/res);
        net.param.s_n2j=round(s_n2j/res);
        net.param.s_dia_x=round(s_dia_x/res);
        net.param.s_dia_y=round(s_dia_y/res);        
        net.param.write_figures=write_figures;
        net.param.write_netlist=write_netlist;
        net.param.fi_req=fi_req;
        net.param.small_edges=0;
        net.param.z_en=z_en;
        net.param.EIT_spacing=1;
        net.param.sim_speed=sim_speed;
        net.param.mode=mode;       
        net.param.ac_freq=ac_freq;
        net.param.debug=debug;
        net.param.spacing_arr=spacing_arr_in/res;
        net.param.dist_step_count=dist_step_count;
        tindex=tindex+1;
        fprintf('T%2g, I1=(%2g,%2g) I2=(%2g,%2g), Zart=(%2g,%2g)',tindex,net.param.i_src_n1i,net.param.i_src_n1j,net.param.i_src_n2i,net.param.i_src_n2j,net.param.r_var_n1k,net.param.r_var_dia);
        
        
        r_dc_arr=imp_image_gen_fn(r_dc,net);
        net.r_dc_arr=r_dc_arr;
        ri_dc_arr=imp_image_gen_fn(ri_dc,net);
        net.ri_dc_arr=ri_dc_arr;
        ci_dc_arr=imp_image_gen_fn(ci_dc,net);
        net.ci_dc_arr=ci_dc_arr;
        % Run
        net=spice_netlist_3d_run_fn(net);
        % Sparse Spice Output
        if net.param.debug==0
            net=read_spice_out_fast_gnd_fn(net);
        end       
        fprintf('\n');
        net_arr2(t0)=net;
    end
    if net.param.debug==0
        net_test=spice_netlist_3d_pp_fn(net_arr2(1));
        if net.param.mode==1
            spice_netlist_3d_plot_ac_fn(net_test);
        else
            spice_netlist_3d_plot_fn(net_test);
        end
        
        if net.param.mode==1 && length(time_arr)>6
            r_var_arr=cat(1,[net_arr2.r_var]);
            r_art_all_arr_arr=[r_var_arr.r_art_all_arr];
            V_node_arr=cat(5,net_arr2.V_node);
            V_node_arr=permute(V_node_arr,[1 2 3 5 4]);
            
            net_t_r=net_arr2(2);
            net_t_r.param.mode=0;
            net_t_i=net_arr2(2);
            net_t_i.param.mode=0;
            %net_t_r.param.sim_speed=3;
            %net_t_i=net_t_r;
            fprintf('Post processing trial %g\n',i')
            for i=1:size(V_node_arr,5)
                % Real
                fprintf('.')
                net_t_r.V_node=real(V_node_arr(:,:,:,:,i));
                net_t_r.t_var=time_arr;
                % Post Processing
                net_t_r=spice_netlist_3d_pp_fn(net_t_r);
                % Plot
                if i==1
                    net_t_r=spice_netlist_3d_plot_fn(net_t_r);
                end
                net_t_r_arr(i)=net_t_r;
                
                % Imaginary
                net_t_i.V_node=imag(V_node_arr(:,:,:,:,i));
                net_t_i.t_var=time_arr;
                % Post Processing
                net_t_i=spice_netlist_3d_pp_fn(net_t_i);              
                net_t_i_arr(i)=net_t_i;
            end
            fprintf('\n')
            net_t_r_t1_arr(t1,:)=net_t_r_arr;
            net_t_i_t1_arr(t1,:)=net_t_i_arr;
        end
    end
end

figure
for i=1:size(net_t_i_t1_arr,1)
    net_t_r_arr=net_t_r_t1_arr(i,:);
    net_t_i_arr=net_t_i_t1_arr(i,:);
    V_s_r_arr=cat(1,net_t_r_arr.V_s);
    V_s_i_arr=cat(1,net_t_i_arr.V_s);
    
    subplot(2,3,1)
    grid on;hold on;    
    plot(net_t_r.f_var(:,1)*1e-3,V_s_r_arr(:,1),'linewidth',2);
    xlabel({'Frequency(kHz)';'(a)'})
    ylabel('|R(V_{DC})| (mV)')
    ax = gca; ax.FontSize=13;
    set(gca, 'FontName', 'Times New Roman');
    legend('1 cm','2 cm','3 cm','fontsize',11)
    subplot(2,3,2)
    grid on;hold on;
    plot(net_t_i.f_var(:,1)*1e-3,abs(V_s_i_arr(:,1)),'linewidth',2);
    xlabel({'Frequency(kHz)';'(b)'})
    ylabel('|I(V_{DC})| (mV)')
    legend('1 cm','2 cm','3 cm','fontsize',11)
    ax = gca; ax.FontSize=13;
    set(gca, 'FontName', 'Times New Roman');
    subplot(2,3,3)
    grid on;hold on;
    plot(V_s_r_arr(:,1),abs(V_s_i_arr(:,1)),'-o','linewidth',2);
    xlabel({'|R(V_{DC})| (mV)';'(c)'})
    ylabel('|I(V_{DC})| (mV)')
    ax = gca; ax.FontSize=13;
    set(gca, 'FontName', 'Times New Roman');
    subplot(2,3,4)
    grid on;hold on;
    plot(net_t_r.f_var(:,1)*1e-3,[net_t_r_arr.V_s_diff],'linewidth',2);
    xlabel({'Frequency(kHz)';'(d)'})
    ylabel('|R(\DeltaV_{PP})| (mV)')
    legend('1 cm','2 cm','3 cm','fontsize',11)
    ax = gca; ax.FontSize=13;
    set(gca, 'FontName', 'Times New Roman');
    subplot(2,3,5)
    grid on;hold on;
    plot(net_t_i.f_var(:,1)*1e-3,abs([net_t_i_arr.V_s_diff]),'linewidth',2);
    xlabel({'Frequency(kHz)';'(e)'})
    ylabel('|I(\DeltaV_{PP})| (mV)')
    legend('1 cm','2 cm','3 cm','fontsize',11)
    ax = gca; ax.FontSize=13;
    set(gca, 'FontName', 'Times New Roman');
    subplot(2,3,6)
    grid on;hold on;
    plot([net_t_r_arr.V_s_diff],[net_t_i_arr.V_s_diff],'-o','linewidth',2);
    xlabel({'|R(\DeltaV_{PP})| (mV)';'(f)'})
    ylabel('|I(\DeltaV_{PP})| (mV)')
    ax = gca; ax.FontSize=13;
    set(gca, 'FontName', 'Times New Roman');
end
set(gcf, 'Position',[159   111   946   539])
file_str=strcat(net_t_r_arr(1).param.path,write_netlist);
print(strcat(file_str,'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
savefig(strcat(file_str,'.fig')); %<-Save as PNG with 300 DP

data=[4.46429	48.79	48.72;
5.20833	48.73	44.95;
6.25	48.45	45.35;
7.8125	48.21	42.08;
8.52273	48.02	43.76;
9.375	48.05	42.84;
10.4167	48.02	41.49;
11.7188	47.18	38.94];


figure
subplot(1,2,1)
grid on;hold on;
subplot1=plot(data(:,1),data(:,2),'o','linewidth',2);
xplot1 = linspace(data(1,1), data(end,1));
% Find coefficients for polynomial (order = 1)
fitResults1 = polyfit(data(:,1),data(:,2),1);
% Evaluate polynomial
yplot1 = polyval(fitResults1,xplot1);
plot(xplot1,yplot1,'--k')
xlabel({'Frequency(kHz)';'(a)'})
ylabel('|R(V_{DC})| (mv)')
set(gca, 'FontName', 'Times New Roman');
ax = gca; ax.FontSize=13;
subplot(1,2,2)
grid on;hold on;
subplot2=plot(data(:,1),data(:,3),'o','linewidth',2);
xplot2 = linspace(data(1,1), data(end,1));
% Find coefficients for polynomial (order = 1)
fitResults2 = polyfit(data(:,1),data(:,3),1);
% Evaluate polynomial
yplot2 = polyval(fitResults2,xplot2);
plot(xplot2,yplot2,'--k')
xlabel({'Frequency(kHz)';'(b)'})
ylabel('|R(\DeltaV_{PP})| (\muv)')
set(gca, 'FontName', 'Times New Roman');
ax = gca; ax.FontSize=13;

set(gcf, 'Position',[403   425   642   241])
file_str=strcat(net_t_r_arr(1).param.path,write_netlist,'_data');
print(strcat(file_str,'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
savefig(strcat(file_str,'.fig')); %<-Save as PNG with 300 DP


figure1=figure;
for i=1:size(net_t_i_t1_arr,1)
    net_t_r_arr=net_t_r_t1_arr(i,:);
    net_t_i_arr=net_t_i_t1_arr(i,:);
    V_s_r_arr=cat(1,net_t_r_arr.V_s);
    V_s_i_arr=cat(1,net_t_i_arr.V_s);
    
    subplot(2,3,2)
    grid on;hold on;    
    plot(net_t_r.f_var(:,1)*1e-3,V_s_r_arr(:,1),'linewidth',2);
    xlabel({'Frequency(kHz)';'(b)'})
    ylabel('|R(V_{DC})| (mV)')
    ax = gca; ax.FontSize=13;
    set(gca, 'FontName', 'Times New Roman');
    legend('1 cm','2 cm','3 cm','fontsize',11)
    subplot(2,3,1)
    grid on;hold on;
    plot(V_s_r_arr(:,1),abs(V_s_i_arr(:,1)),'-o','linewidth',2);
    xlabel({'|R(V_{DC})| (mV)';'(a)'})
    ylabel('|I(V_{DC})| (mV)')
    ax = gca; ax.FontSize=13;
    set(gca, 'FontName', 'Times New Roman');
    subplot(2,3,5)
    grid on;hold on;
    plot(net_t_r.f_var(:,1)*1e-3,[net_t_r_arr.V_s_diff],'linewidth',2);
    xlabel({'Frequency(kHz)';'(e)'})
    ylabel('|R(\DeltaV_{PP})| (mV)')
    legend('1 cm','2 cm','3 cm','fontsize',11)
    ax = gca; ax.FontSize=13;
    set(gca, 'FontName', 'Times New Roman');
    subplot(2,3,4)
    grid on;hold on;
    plot([net_t_r_arr.V_s_diff],[net_t_i_arr.V_s_diff],'-o','linewidth',2);
    xlabel({'|R(\DeltaV_{PP})| (mV)';'(d)'})
    ylabel('|I(\DeltaV_{PP})| (mV)')
    ax = gca; ax.FontSize=13;
    set(gca, 'FontName', 'Times New Roman');
end

axes1 = axes('Parent',figure1,...
    'Position',[0.714876544860762 0.645925347650311 0.22677250376503 0.279074652349689]);
hold(axes1,'on');

% Create plot
plot(data(:,1),data(:,2),'Parent',axes1,'Marker','o','LineWidth',2,'LineStyle','none');

% Create plot
plot(xplot1,yplot1,'Parent',axes1,'LineStyle','--','Color',[0 0 0]);

% Create ylabel
ylabel('|R(V_{DC})| (mv)','FontName','Times New Roman');

% Create xlabel
xlabel({'Frequency(kHz)','(c)'},'FontName','Times New Roman');

grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',13);
% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.714876544860762 0.172088138347986 0.233537831460591 0.279074652349689]);
hold(axes2,'on');

% Create plot
plot(data(:,1),data(:,3),'Parent',axes2,'Marker','o','LineWidth',2,'LineStyle','none');

% Create plot
plot(xplot2,yplot2,'Parent',axes2,'LineStyle','--','Color',[0 0 0]);

% Create ylabel
ylabel('|R(\DeltaV_{PP})| (\muv)','FontName','Times New Roman');

% Create xlabel
xlabel({'Frequency(kHz)','(f)'},'FontName','Times New Roman');

grid(axes2,'on');
% Set the remaining axes properties
set(axes2,'FontName','Times New Roman','FontSize',13);

set(gcf, 'Position',[159   111   946   539])
file_str=strcat(net_t_r_arr(1).param.path,write_netlist);
print(strcat(file_str,'_all.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
savefig(strcat(file_str,'_all.fig')); %<-Save as PNG with 300 DP

if write_mat==1 && debug==0
    if ~exist(net.param.path,'dir')
        mkdir(net.param.path);
    end
    save(strcat(net.param.path,'netlist.mat'));
end

