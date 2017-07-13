function out=obstacle(t)
global obs1_path 
t_data=(0:200)'/10;
out=interp1(t_data,[obs1_path(1,:) ;obs1_path],t);

end
