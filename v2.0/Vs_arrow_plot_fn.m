function Vs_arrow_plot_fn(net)

for i=1:length(net.param.s_n1i)
    %plot([net.param.s_n1j(i)+1,net.param.s_n2j(i)+1],[net.param.s_n1i(i)+1,net.param.s_n2i(i)+1],'linewidth',2)
    %text(net.param.s_n2j(i)+1,net.param.s_n2i(i)+1,'\leftarrow t_1','FontSize',12,'FontWeight','bold')
    %annotation('textarrow',[net.param.s_n1j(i)+1,net.param.s_n2j(i)+1]./(net.Nj+1),[net.param.s_n1i(i)+1,net.param.s_n2i(i)+1]./(net.Ni+1),'String','y = x ')
    
    drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );  

    x1 = fliplr([net.param.s_n1j(i)+1,net.param.s_n2j(i)+1]);
    y1 = fliplr([net.param.s_n1i(i)+1,net.param.s_n2i(i)+1]);

    drawArrow(x1,y1,'linewidth',2,'color',net.param.color_arr(1,:)); hold on
    %,'color',net.param.color(1,:)
end