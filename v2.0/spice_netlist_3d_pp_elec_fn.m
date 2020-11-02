function net=spice_netlist_3d_pp_elec_fn(net)

%if net(1).param.re_var_max(1)>0
    % Voltage Sensing V elec
    if ~isempty(net.s)
        %if net.param.mode~=1
        for i=1:length(net.s)
            [net.s(i).PAT_e,net.s(i).dVpp_e,net.s(i).Vdc_e]=find_delay_fn(net.t_var,abs(net.s(i).Ve),net.param.r_var_freq(1),0);
            [net.s(i).PAT_e_r,net.s(i).dVpp_e_r,net.s(i).Vdc_e_r]=find_delay_fn(net.t_var,real(net.s(i).Ve),net.param.r_var_freq(1),0);
            [net.s(i).PAT_e_i,net.s(i).dVpp_e_i,net.s(i).Vdc_e_i]=find_delay_fn(net.t_var,imag(net.s(i).Ve),net.param.r_var_freq(1),0);
        end
        %end
    end
    % Current Injection V elec
    if ~isempty(net.curr)
        %if net.param.mode~=1
        for i=1:length(net.curr)
            [net.curr(i).PAT_e,net.curr(i).dVpp_e,net.curr(i).Vdc_e]=find_delay_fn(net.t_var,abs(net.curr(i).Ve),net.param.r_var_freq(1),0);
            [net.curr(i).PAT_e_r,net.curr(i).dVpp_e_r,net.curr(i).Vdc_e_r]=find_delay_fn(net.t_var,real(net.curr(i).Ve),net.param.r_var_freq(1),0);
            [net.curr(i).PAT_e_i,net.curr(i).dVpp_e_i,net.curr(i).Vdc_e_i]=find_delay_fn(net.t_var,imag(net.curr(i).Ve),net.param.r_var_freq(1),0);
        end
        %end
    end
    % Voltage Sensing V skin
    if ~isempty(net.s)
        for i=1:length(net.s)
            for d=1:length(net.s(i).n1i)
                net.s(i).Vs_n1(d,:)=squeeze((net.V_node(net.s(i).n1i(d),net.s(i).n1j(d),net.s(i).n1k(d),:)));
                net.s(i).Vs_n2(d,:)=squeeze((net.V_node(net.s(i).n2i(d),net.s(i).n2j(d),net.s(i).n2k(d),:)));
            end
            net.s(i).Vs_a=net.s(i).Vs_n2-net.s(i).Vs_n1;
            net.s(i).Vs=squeeze(mean(net.s(i).Vs_a,1));
        end
        %if net.param.mode~=1
        for i=1:length(net.s)
            [net.s(i).PAT_s,net.s(i).dVpp_s,net.s(i).Vdc_s]=find_delay_fn(net.t_var,abs(net.s(i).Vs),net.param.r_var_freq(1),0);
            [net.s(i).PAT_s_r,net.s(i).dVpp_s_r,net.s(i).Vdc_s_r]=find_delay_fn(net.t_var,real(net.s(i).Vs),net.param.r_var_freq(1),0);
            [net.s(i).PAT_s_i,net.s(i).dVpp_s_i,net.s(i).Vdc_s_i]=find_delay_fn(net.t_var,imag(net.s(i).Vs),net.param.r_var_freq(1),0);
        end
        %end
    end
    % Current Injection V skin
    if ~isempty(net.curr)
        for i=1:length(net.curr)
            for d=1:length(net.curr(i).n1i)
                net.curr(i).Vs_n1(d,:)=squeeze((net.V_node(net.curr(i).n1i(d),net.curr(i).n1j(d),net.curr(i).n1k(d),:)));
                net.curr(i).Vs_n2(d,:)=squeeze((net.V_node(net.curr(i).n2i(d),net.curr(i).n2j(d),net.curr(i).n2k(d),:)));
            end
            net.curr(i).Vs_a=net.curr(i).Vs_n2-net.curr(i).Vs_n1;
            net.curr(i).Vs=squeeze(mean(net.curr(i).Vs_a,1));
        end
        %if net.param.mode~=1
        for i=1:length(net.curr)
            [net.curr(i).PAT_s,net.curr(i).dVpp_s,net.curr(i).Vdc_s]=find_delay_fn(net.t_var,abs(net.curr(i).Vs),net.param.r_var_freq(1),0);
            [net.curr(i).PAT_s_r,net.curr(i).dVpp_s_r,net.curr(i).Vdc_s_r]=find_delay_fn(net.t_var,real(net.curr(i).Vs),net.param.r_var_freq(1),0);
            [net.curr(i).PAT_s_i,net.curr(i).dVpp_s_i,net.curr(i).Vdc_s_i]=find_delay_fn(net.t_var,imag(net.curr(i).Vs),net.param.r_var_freq(1),0);
        end
        %end
    end
%end