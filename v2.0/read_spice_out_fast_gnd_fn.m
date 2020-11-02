function net=read_spice_out_fast_gnd_fn(net)

%   Bassem Ibrahim
%   2020.11.02
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari
%
%   function net=spice_netlist_3d_plot_ac_fn(net)
%   Read and sparse SPICE output file
%   net.V_node: Node voltage
%   net.f_var: Frequency array
%   net.t_var: Time array
tic
spice_out_file=net(1).param.output_file;
header = importfile_ltspice_header(spice_out_file);
str=table2array(header(3,:));
[a0 b0]=regexp(str,'DC');
[a1 b1]=regexp(str,'Transient');
[a2 b2]=regexp(str,'AC');
% DC Analysis
if ~isempty(a0)
    test = importfile_ltspice_raw(spice_out_file);
    step_index=find(~isnan(table2array(test(:,1))));
    voltage_index_raw=regexp(test.VarName3,'V(\d+');
    voltage_index=[];
    for i=1:length(voltage_index_raw)
        if ~isempty(voltage_index_raw{i})
            voltage_index=[voltage_index;i];
        end
    end
    voltage_field = strrep(test.VarName3(voltage_index),'(','');
    voltage_field = strrep(voltage_field,')','');
    voltage_node_index = str2num(char(strrep(voltage_field,'V','')));
    [a b]=sort(voltage_node_index);
    L=length(step_index);
    for i=1:L
        V_raw=[test.VarName2(step_index(i)+1:step_index(i)+length(voltage_index))];
        V=[V_raw(b)];
        V_node(:,:,:,i)=permute(reshape(V,net.Ni+1,net.Nj+1,net.Nk+1),[2 1 3]);
        t_var(i,1)=str2num(char(test.VarName3(step_index(i))));
    end
    net.V_node=V_node;
    net.t_var=t_var;
    
    % Transient Analysis
elseif ~isempty(a1)
    test = importfile_ltspice_raw_fast_fn(spice_out_file);
    step_index=find(table2array(test(:,1))~="");
    voltage_index_raw=regexp(test.VarName3,'V\(\d+\)');
    current_index_raw=regexp(test.VarName3,'I\(I\d+\)');
    voltage_index=[];
    current_index=[];
    for i=1:length(voltage_index_raw)
        if ~isempty(voltage_index_raw{i})
            voltage_index=[voltage_index;i];
        end
        if ~isempty(current_index_raw{i})
            current_index=[current_index;i];
        end
    end
    voltage_field = strrep(test.VarName3(voltage_index),'(','');
    voltage_field = strrep(voltage_field,')','');
    voltage_node_index = str2num(char(strrep(voltage_field,'V','')));
    [a b]=sort(voltage_node_index);
    step_index=step_index(2:end);
    L=length(step_index);
    for i=1:L
        V_raw=[test.VarName2(step_index(i)-1+voltage_index)];
        V=[V_raw(b)];
        V_node(:,:,:,i)=permute(reshape(V,net.Nj+1,net.Ni+1,net.Nk+1),[2 1 3]);
        t_var(i,1)=str2num(char(test.VarName3(step_index(i))));
    end
    net.V_node=V_node;
    net.t_var=t_var;
    net.f_var=0;
    % AC Analysis
elseif ~isempty(a2)
    test = importfile_ltspice_raw_fast_fn(spice_out_file);
    step_index=find(table2array(test(:,1))~="");
    voltage_index_raw=regexp(test.VarName3,'V\(\d+\)');
    voltage_index=[];
    for i=1:length(voltage_index_raw)
        if ~isempty(voltage_index_raw{i})
            voltage_index=[voltage_index;i];
        end
    end
    voltage_field = strrep(test.VarName3(voltage_index),'(','');
    voltage_field = strrep(voltage_field,')','');
    voltage_node_index = str2num(char(strrep(voltage_field,'V','')));
    [a b]=sort(voltage_node_index);
    step_index=step_index(2:end);
    L=length(step_index);
    for i=1:L
        Vr_raw=[test.VarName2(step_index(i)-1+voltage_index)];
        newStr = split(test.VarName3(step_index(i)-1+voltage_index),',');
        Vi_raw=[str2num(char(newStr{:,2}))];
        V_raw=Vr_raw+1i*Vi_raw;
        V=[V_raw(b)];
        V_node=permute(reshape(V,net(1).Nj+1,net(1).Ni+1,net(1).Nk+1),[2 1 3]);
        net(i).V_node_x=V_node;
        for m=1:length(net(1).curr)
            [r c]=find(voltage_node_index==net(1).curr(m).n1);
            net(i).curr(m).Ve_n1_x=V_raw(r);
            [r c]=find(voltage_node_index==net(1).curr(m).n2);
            net(i).curr(m).Ve_n2_x=V_raw(r);
            net(i).curr(m).Ve_a_x=net(i).curr(m).Ve_n2_x-net(i).curr(m).Ve_n1_x;
            net(i).curr(m).Ve_x=net(i).curr(m).Ve_a_x(1,1);
        end
        %net(20).curr(1).Vs_n2(1,:)=squeeze((net(20).V_node_x(net(1).curr(1).n2i(1),net(1).curr(1).n2j(1),1)));
        for m=1:length(net(1).s)
            [r c]=find(voltage_node_index==net(1).s(m).n1);
            net(i).s(m).Ve_n1_x=V_raw(r);
            [r c]=find(voltage_node_index==net(1).s(m).n2);
            net(i).s(m).Ve_n2_x=V_raw(r);
            net(i).s(m).Ve_a_x=net(i).s(m).Ve_n2_x-net(i).s(m).Ve_n1_x;
            net(i).s(m).Ve_x=net(i).s(m).Ve_a_x(1,1);
        end
        %net(20).s(1).Vs_n2(1,:)=squeeze((net(20).V_node_x(net(1).s(1).n2i(1),net(1).s(1).n2j(1),1)));
        newStr = split(test.VarName3(step_index(i)),',');
        f_var(i,1)=str2num(char(newStr{1}));
        f_var(i,2)=str2num(char(newStr{2}));
    end
    net(1).f_var=f_var;
    net(1).t_var=0;
    
    
    % Get electrodes voltage
    if net(1).param.se_en==1
        voltage_se_index_raw=regexp(test.VarName3,'V\(\d+\_se)');
        voltage_se_index=[];
        for i=1:length(voltage_se_index_raw)
            if ~isempty(voltage_se_index_raw{i})
                voltage_se_index=[voltage_se_index;i];
            end
        end
        voltage_se_field = strrep(test.VarName3(voltage_se_index),'(','');
        voltage_se_field = strrep(voltage_se_field,')','');
        voltage_se_node_index = ((strrep(voltage_se_field,'V','')));
        voltage_se_node_index = str2num(char(strrep(voltage_se_node_index,'_se','')));
        for i=1:L
            Vr_se_raw=[test.VarName2(step_index(i)-1+voltage_se_index)];
            newStr_se = split(test.VarName3(step_index(i)-1+voltage_se_index),',');
            Vi_se_raw=[str2num(char(newStr_se{:,2}))];
            V_se_raw=Vr_se_raw+1i*Vi_se_raw;
            for m=1:length(net(1).curr)
                [r c]=find(voltage_se_node_index==net(1).curr(m).n1);
                net(i).curr(m).Ve_n1_x=V_se_raw(r,:);
                [r c]=find(voltage_se_node_index==net(1).curr(m).n2);
                net(i).curr(m).Ve_n2_x=V_se_raw(r,:);
                net(i).curr(m).Ve_a_x=net(i).curr(m).Ve_n2_x-net(i).curr(m).Ve_n1_x;
                net(i).curr(m).Ve_x=net(i).curr(m).Ve_a_x(1,1);
            end
            for m=1:length(net(1).s)
                [r c]=find(voltage_se_node_index==net(1).s(m).n1);
                net(i).s(m).Ve_n1_x=V_se_raw(r,:);
                [r c]=find(voltage_se_node_index==net(1).s(m).n2);
                net(i).s(m).Ve_n2_x=V_se_raw(r,:);
                net(i).s(m).Ve_a_x=net(i).s(m).Ve_n2_x-net(i).s(m).Ve_n1_x;
                net(i).s(m).Ve_x=net(i).s(m).Ve_a_x(1,1);
            end
        end
    end
    
    % current
    if net(1).param.icalc_en==1
        if net(1).param.z_en==0
            current_index_raw=regexp(test.VarName3,'I\(.*_\d+_\d+\)');
            current_index=[];
            for i=1:length(current_index_raw)
                if ~isempty(current_index_raw{i})
                    current_index=[current_index;i];
                end
            end
            
            for i=1:length(current_index)
                str=test.VarName3(current_index(i));
                [a b] = regexpi(str,'_(\d+)_(\d+)','match','tokens');
                current_field1(i)=str2num(b{1}{1});
                current_field2(i)=str2num(b{1}{2});
            end
            
            for i=1:L
                Ir_raw=[test.VarName2(step_index(i)-1+current_index)];
                newStr = split(test.VarName3(step_index(i)-1+current_index),',');
                Ii_raw=[str2num(char(newStr{:,2}))];
                I_raw=Ir_raw+1i*Ii_raw;
                
                I_node.x.Ri=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.y.Ri=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.z.Ri=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.x.Re=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.y.Re=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.z.Re=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.x.Ci=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.y.Ci=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.z.Ci=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                
                for n=1:length(I_raw)
                    [node_idx(1),node_idx(2),node_idx(3)]=ind2sub([net(1).Ni+1,net(1).Nj+1,net(1).Nk+1],find(current_field1(n)==net(1).node_arr_raw));
                    if node_idx(1)==net(1).Ni+1
                        if node_idx(2)==net(1).Nj+1 && current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2),node_idx(3)+1)
                            I_node.z.Re(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                        else
                            if current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2)+1,node_idx(3))
                                I_node.x.Re(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                            elseif current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2),node_idx(3)+1)
                                I_node.z.Re(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                            else
                                fprintf('Error in current nodes\n');
                                break;
                            end
                        end
                    elseif node_idx(2)==net(1).Nj+1
                        if node_idx(3)==net(1).Nk+1 && current_field2(n)==net(1).node_arr_raw(node_idx(1)+1,node_idx(2),node_idx(3))
                            I_node.y.Re(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                        else
                            if current_field2(n)==net(1).node_arr_raw(node_idx(1)+1,node_idx(2),node_idx(3))
                                I_node.y.Re(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                            elseif current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2),node_idx(3)+1)
                                I_node.z.Re(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                            else
                                fprintf('Error in current nodes\n');
                                break;
                            end
                        end
                    elseif node_idx(3)==net(1).Nk+1
                        if node_idx(1)==net(1).Ni+1 && current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2)+1,node_idx(3))
                            I_node.x.Re(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                        else
                            if current_field2(n)==net(1).node_arr_raw(node_idx(1)+1,node_idx(2),node_idx(3))
                                I_node.y.Re(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                            elseif current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2)+1,node_idx(3))
                                I_node.x.Re(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                            else
                                fprintf('Error in current nodes\n');
                                break;
                            end
                        end
                    else
                        if current_field2(n)==net(1).node_arr_raw(node_idx(1)+1,node_idx(2),node_idx(3))
                            I_node.y.Re(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                        elseif current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2)+1,node_idx(3))
                            I_node.x.Re(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                        elseif current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2),node_idx(3)+1)
                            I_node.z.Re(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                        else
                            fprintf('Error in current nodes\n');
                            break;
                        end
                    end
                end
                net(i).I_node=I_node;
                %         I_node_arr(:,:,:,1)=I_node.x;
                %         I_node_arr(:,:,:,2)=I_node.y;
                %         I_node_arr(:,:,:,3)=I_node.z;
            end
        else
            current_index_raw=regexp(test.VarName3,'I\(.*_\d+_\d+.*\)');
            current_index=[];
            for i=1:length(current_index_raw)
                if ~isempty(current_index_raw{i})
                    current_index=[current_index;i];
                end
            end
            
            for i=1:length(current_index)
                str=test.VarName3(current_index(i));
                [a b] = regexpi(str,'(\w)_(\d+)_(\d+)_(\w)','match','tokens');
                current_elem{i}=(b{1}{1});
                current_field1(i)=str2num(b{1}{2});
                current_field2(i)=str2num(b{1}{3});
                current_loc{i}=(b{1}{4});
            end
            
            for i=1:L
                Ir_raw=[test.VarName2(step_index(i)-1+current_index)];
                newStr = split(test.VarName3(step_index(i)-1+current_index),',');
                Ii_raw=[str2num(char(newStr{:,2}))];
                I_raw=Ir_raw+1i*Ii_raw;
                
                I_node.x.Ri=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.y.Ri=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.z.Ri=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.x.Re=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.y.Re=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.z.Re=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.x.Ci=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.y.Ci=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                I_node.z.Ci=zeros(net(1).Ni+1,net(1).Nj+1,net(1).Nk+1);
                for n=1:length(I_raw)
                    [node_idx(1),node_idx(2),node_idx(3)]=ind2sub([net(1).Ni+1,net(1).Nj+1,net(1).Nk+1],find(current_field1(n)==net(1).node_arr_raw));
                    if node_idx(1)==net(1).Ni+1
                        if node_idx(2)==net(1).Nj+1 && current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2),node_idx(3)+1)
                            I_node.z.(strcat(current_elem{n},current_loc{n}))(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                        else
                            if current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2)+1,node_idx(3))
                                I_node.x.(strcat(current_elem{n},current_loc{n}))(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                            elseif current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2),node_idx(3)+1)
                                I_node.z.(strcat(current_elem{n},current_loc{n}))(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                            else
                                fprintf('Error in current nodes\n');
                                break;
                            end
                        end
                    elseif node_idx(2)==net(1).Nj+1
                        if node_idx(3)==net(1).Nk+1 && current_field2(n)==net(1).node_arr_raw(node_idx(1)+1,node_idx(2),node_idx(3))
                            I_node.y.(strcat(current_elem{n},current_loc{n}))(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                        else
                            if current_field2(n)==net(1).node_arr_raw(node_idx(1)+1,node_idx(2),node_idx(3))
                                I_node.y.(strcat(current_elem{n},current_loc{n}))(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                            elseif current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2),node_idx(3)+1)
                                I_node.z.(strcat(current_elem{n},current_loc{n}))(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                            else
                                fprintf('Error in current nodes\n');
                                break;
                            end
                        end
                    elseif node_idx(3)==net(1).Nk+1
                        if node_idx(1)==net(1).Ni+1 && current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2)+1,node_idx(3))
                            I_node.x.(strcat(current_elem{n},current_loc{n}))(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                        else
                            if current_field2(n)==net(1).node_arr_raw(node_idx(1)+1,node_idx(2),node_idx(3))
                                I_node.y.(strcat(current_elem{n},current_loc{n}))(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                            elseif current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2)+1,node_idx(3))
                                I_node.x.(strcat(current_elem{n},current_loc{n}))(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                            else
                                fprintf('Error in current nodes\n');
                                break;
                            end
                        end
                    else
                        if current_field2(n)==net(1).node_arr_raw(node_idx(1)+1,node_idx(2),node_idx(3))
                            I_node.y.(strcat(current_elem{n},current_loc{n}))(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                        elseif current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2)+1,node_idx(3))
                            I_node.x.(strcat(current_elem{n},current_loc{n}))(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                        elseif current_field2(n)==net(1).node_arr_raw(node_idx(1),node_idx(2),node_idx(3)+1)
                            I_node.z.(strcat(current_elem{n},current_loc{n}))(node_idx(1),node_idx(2),node_idx(3))=I_raw(n);
                        else
                            fprintf('Error in current nodes\n');
                            break;
                        end
                    end
                end
                net(i).I_node=I_node;
                %         I_node_arr(:,:,:,1)=I_node.x;
                %         I_node_arr(:,:,:,2)=I_node.y;
                %         I_node_arr(:,:,:,3)=I_node.z;
            end
        end
    else
        net(i).I_node=0;
    end
end
fprintf('Output files sparsed:');toc;