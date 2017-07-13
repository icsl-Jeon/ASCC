function out=xdot_ee(t)
global xdotl_repro

t_data=(0:200)'/10;
out=interp1(t_data,[xdotl_repro(:,1) xdotl_repro]'/10,t);

end
