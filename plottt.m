% plotting 
q_sim=interp1(output1(:,end)/t_scale,next_q.Data(:,[1 3 4 5]),1:nbDatab);
figure()

for i=1:3
    subplot(3,1,i)
    xlabel('number of Data')
    tit=strcat('q',num2str(i));
    ylabel(tit)
    plot(1:200,q_sim(:,i+1))
    hold on
    plot(1:3:200,qs(1:3:200,i+2),'k.')
    hold off
    legend('acheieved trajectory','learned joint trajectory')
    gtext(tit)
end
gtext('number of data')
