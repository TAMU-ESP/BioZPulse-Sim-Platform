function net=read_spice_out_fast_gnd_fn(net)

%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling,
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   function net=spice_netlist_3d_plot_ac_fn(net)
%   Read and sparse SPICE output file
%   net.V_node: Node voltage
%   net.f_var: Frequency array
%   net.t_var: Time array

spice_out_file=net.param.output_file;
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
        V_raw=[test.VarName2(step_index(i)+1:step_index(i)+length(voltage_index))]*1e3;
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
        V_raw=[test.VarName2(step_index(i)-1+voltage_index)]*1e3;
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
        Vr_raw=[test.VarName2(step_index(i)-1+voltage_index)]*1e3;
        newStr = split(test.VarName3(step_index(i)-1+voltage_index),',');
        Vi_raw=[str2num(char(newStr{:,2}))]*1e3;
        V_raw=Vr_raw+1i*Vi_raw;
        V=[V_raw(b)];
        V_node(:,:,:,i)=permute(reshape(V,net.Nj+1,net.Ni+1,net.Nk+1),[2 1 3]);               
        newStr = split(test.VarName3(step_index(i)),',');
        f_var(i,1)=str2num(char(newStr{1}));
        f_var(i,2)=str2num(char(newStr{2}));
    end
    
    net.V_node=V_node;
    net.f_var=f_var;
    net.t_var=0;
end
