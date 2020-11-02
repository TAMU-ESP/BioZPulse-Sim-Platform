function txt = myupdatefcn(~,event_obj, A)
% Customizes text of data tips
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
txt = {['X: ',num2str(pos(1))],...
    ['Y: ',num2str(pos(2))],...
    ['U: ',num2str(A(pos(2),pos(1)))]};