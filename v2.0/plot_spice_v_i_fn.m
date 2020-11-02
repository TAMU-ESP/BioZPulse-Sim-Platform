function net=plot_spice_v_i_fn(net,opt,y,V_node_min,V_node_max)

%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling,
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   function net=spice_netlist_3d_plot_fn(net)
%   Plot Output Results


s_arr=[net.s];
Vs_arr=abs([s_arr.Vs]);
Vs_max=max(Vs_arr,[],'all');
figure
subplot(1,2,1)
if Vs_max<10
    mesh_plot_fn(squeeze(abs(net.V_node(y,:,:)))'*1e3,1,[V_node_min V_node_max]*1e3);
    title(strcat('2D Node and Electrode Voltage(Ve) (uV) at Y=',num2str(y)));
else
    mesh_plot_fn(squeeze(abs(net.V_node(y,:,:)))',1,[V_node_min V_node_max]);
    title(strcat('2D Node and Electrode Voltage(Ve) (mV) at Y=',num2str(y)));
end
plot_velec_fn(net,opt);
xlabel('X');ylabel('Z')
ax = gca; ax.FontSize=13;

subplot(1,2,2)
plot_core_3d_fn(A,net,0,[V_node_min V_node_max]*1e3);
if Vs_max<10
    A=abs(net.V_node)*1e3;
    title('3D Node Voltage Map (mV)');
else
    A=abs(net.V_node);
    title('3D Node Voltage Map (uV)');
end
xlabel('Y');ylabel('X');zlabel('Z');
ax = gca; ax.FontSize=13;


