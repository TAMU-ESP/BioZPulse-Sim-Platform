function plot_grid_fn(A)

hold on;
for i = 1:size(A,1)  
   plot([i-.5,i-.5],[.5,size(A,1)+.5],'k-');
end
for i = 1:size(A,2)
   plot([.5,size(A,2)+.5],[i-.5,i-.5],'k-'); 
end