function out=qdot_ee(t)
global qdot_DMP
t_data=(0:200)'/10;
out=interp1(t_data,[qdot_DMP(:,1) qdot_DMP]'/10,t);

end
