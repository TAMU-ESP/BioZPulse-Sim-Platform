%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling,
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   Simulate DC Bio-Z(V_DC) for different horizontal distance of sensing electrodes


clc
clear
close all
fclose('all');
warning('off','all');


debug=0;
write_netlist='dc_X';

% Dimensions
Nk=7;
% Current Source
date_start=datestr(now,'yyyy-mm-dd');
i_num=1;
i_opp=0;
i_freq=100;
dia_ar=1;

% Configurations
z_en=0;
nstep=8;
mode=1;  %0=> .TRAN, 1=> .AC
ac_freq=[1 1e3 100e3]; %npts fstart fend
fi_req=[1,3,4];
sim_speed=4;
dist_step_count=10;
%dia_arr=[6,8];
s_num=21;
dia_arr=[0,4];
L_arr=[2];
for t2=1:length(L_arr)
    for t1=1:length(dia_arr)
        clear net_arr2 net;
        margin_x=max(dia_arr(t1)*3,8);
        margin_y=max(dia_arr(t1)*2,6);
        spacing_x=max(dia_arr(t1)*2,6);
        spacing_y=max(dia_arr(t1)*2,6);
        L=L_arr(t2);        
        % Sensing
        dia_x=dia_arr(t1)*dia_ar;
        dia_y=dia_arr(t1);        
        spacing_arr_in=spacing_y;
        Ni=(max(dia_arr)+(3)*spacing_y+2*margin_y);
        Nj=(max(dia_arr)*dia_ar+(s_num*1-1)*spacing_x+2*margin_x);
        i_dia_x=dia_x;
        s_dia_x=dia_x;
        i_dia_y=dia_y;
        s_dia_y=dia_y;
        i_src_n1i_arr=round((Ni-(dia_y+(1*1+2)*spacing_y))/2);
        i_src_n2i_arr=i_src_n1i_arr+spacing_y*3;
        i_src_n1j_arr=round(Nj/2-dia_x/2);
        i_src_n2j_arr=i_src_n1j_arr;
        s_n1i=[];s_n1j=[];s_n2i=[];s_n2j=[];
        if s_num>0
            for i=1:s_num
                s_n1j(i)=i_src_n1j_arr-(s_num-1)/2*spacing_x+spacing_x*(1*i-1);
                s_n2j(i)=s_n1j(i);
                s_n1i(i)=i_src_n1i_arr+spacing_y;
                s_n2i(i)=i_src_n1i_arr+spacing_y*2;
            end
        end        
        % Artery Model
        r_dc=30;
        ri_dc=r_dc;
        ci_dc=100e-9;
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
        ci_var_max=ci_art_dc*1e-2*ones(1,S);
        r_var_steps=3*ones(1,S);
        r_var_v=0*ones(1,S);       % 1-> vertical, 0-> horizontal
        r_var_dim=3*ones(1,S);     % Dimension (1,2,3)
        r_var_nan=1*ones(1,S);     % 0-> No NaN, 1-> NaN
        r_var_dia=2*ones(1,S);
        
        %time_en=1; %1=> time in R,C var
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
            i_src_n1j=i_src_n1j_arr(1);
            r_var_delay=r_var_delay_arr(1);
            r_var_freq=r_var_freq_arr(1,:);
            r_var_time=r_var_time_arr(1);
            tstep=r_var_time*nstep;
            i_src_n1i=i_src_n1i_arr(1);
            i_src_n2i=i_src_n2i_arr(1);
            r_var_n1k=r_var_n1k_arr(1);
            time=time_arr(t0);
            
            net.Ni=round(Ni/L);
            net.Nj=round(Nj/L);
            net.Nk=round(Nk/L);
            net.param.i_src_n1i=round(i_src_n1i/L);
            net.param.i_src_n2i=round(i_src_n2i/L);
            net.param.i_src_n1j=round(i_src_n1j/L);
            net.param.i_src_n2j=round(i_src_n1j/L);
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
            net.param.r_var_n1i=round(r_var_n1i/L);
            net.param.r_var_n1j=round(r_var_n1j/L);
            net.param.r_var_n1k=round(r_var_n1k/L);
            net.param.r_var_v=r_var_v;
            net.param.r_var_dim=r_var_dim;
            net.param.r_var_num=net.Ni*ones(1,S);
            net.param.r_var_time=r_var_time;
            net.param.r_var_seg=net.Ni*ones(1,S);
            net.param.r_var_dia=round(r_var_dia/L);
            net.param.r_var_delay=r_var_delay;
            net.param.r_var_freq=r_var_freq;
            net.param.tstep=tstep;
            net.param.time=time;
            net.curr(1).dia_x=round(i_dia_x/L);
            net.curr(1).dia_y=round(i_dia_y/L);
            net.curr(1).ival=i_val;
            net.curr(1).iac=i_ac;
            net.curr(1).freq=i_freq;
            net.param.s_n1i=round(s_n1i/L);
            net.param.s_n2i=round(s_n2i/L);
            net.param.s_n1j=round(s_n1j/L);
            net.param.s_n2j=round(s_n2j/L);
            net.param.s_dia_x=round(s_dia_x/L);
            net.param.s_dia_y=round(s_dia_y/L);
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
            net.param.spacing_arr=spacing_arr_in/L;
            net.param.dist_step_count=dist_step_count;
            net.param.L=L;
            net.param.date=date_start;
            
            tindex=tindex+1;
            fprintf('Time step %2g\n',tindex);
            
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
                net_t_r.param.sim_speed=3;
                
                for i=1:size(V_node_arr,5)
                    % Real
                    fprintf('Post processing trial %g\n',i')
                    net_t_r.V_node=real(V_node_arr(:,:,:,:,i));
                    net_t_r.t_var=time_arr;
                    % Post Processing
                    net_t_r=spice_netlist_3d_pp_fn(net_t_r);
                    % Plot
                    if i==1
                        net_t_r=spice_netlist_3d_plot_fn(net_t_r);
                    end
                    net_t_r_arr(i)=net_t_r;
                end
                net_t_i_t1_arr(t2,t1,:)=net_t_r_arr;
            end
        end
    end
end



