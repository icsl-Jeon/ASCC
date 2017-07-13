%% parameters 

% robot parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global arms dof arms_flat

% potential field parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global alpha0_max eta po beta lambda 
dof=5;



% DMP learning parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbDatal=200; nbVarl=dof; alphal=1;  nbStatesl=10; Kl=500; Dl=40; 
dt=0.01;
nbDatab=200; nbVarb=3; alphab=1;  nbStatesb=12; Kb=200; Db=20; 



%% Link defined

LinkLength=[1 1 4 4 5];

arms=cell(dof,1); arms_flat=cell(2*dof,1);
eta=800; % need to be tuned
alpha0_max=3;
po=3;
beta=3;
lambda=40;

[arms,arms_flat]=Make_arm(LinkLength,SE3(eye(4)));
[arms0,arms_flat0]=Make_arm(LinkLength,SE3(eye(4)));
%% Body path & obstacle defined 
% straiht line 
% 
xb(1,:)=linspace(25,-25,nbDatab);
xb(2,:)=zeros(1,nbDatab);
xb(3,:)=5*ones(1,nbDatab);

% obstacle path generation  

xo10=[15 10 -4];
xo1g=[15 -10 -4];


waypoint1=xb(:,nbDatab/2)'+[0 0 -3];


obs_path1=mtraj(@lspb,xo10,waypoint1,nbDatab/4);
obs_path2=mtraj(@lspb,waypoint1,xo1g,nbDatab/4*3);




% % C -curve
% 
% th=linspace(-pi,0,nbDatab);
% for i=1:nbDatab
% xb(1,i)=25*cos(th(i));
% xb(2,i)=25*sin(th(i));
% end
% 
% xb(3,:)=5*ones(1,nbDatab);
% 
% % ¿Ï¸¸ÇÑ °î¼± 
% xb(2,:)=xb(2,:)/2;



global obs1_path obs2

% obs1_path=repmat(xb(:,nbDatab/2)',nbDatab,1)+[0 0 -3];
% obs1_path=[obs_path1 ; obs_path2];
obs1_path=[repmat(waypoint1(1),200,1)  [linspace(-10,0,100)  linspace(0.03,10,100)]'  repmat(waypoint1(1),200,1)]; 


obs2=xb(:,nbDatab*3/4)+[0 100 100]';





xb0=xb(:,1);
xbg=xb(:,end);

[wb,Mub,Sigmab]=DMP_learning(nbDatab,nbVarb,alphab,dt,nbStatesb,Kb,Db,xb);
[xb_repro,xdotb_repro]=DMP_repro(nbDatab,3,alphab,dt,nbStatesb,Kb,Db,wb,Mub,Sigmab,1,xb0,xb0,xbg,0);


global yaw
% yaw calculation 
yaw=zeros(1,nbDatab);
yaw(1:end-1)=atan2(xdotb_repro(2,2:end),xdotb_repro(1,2:end)); yaw(end)=yaw(end-1);
yaw(yaw<0)=yaw(yaw<0)+2*pi;
for check=1:nbDatab-1
    if abs(yaw(check+1)-yaw(check))>pi
        yaw(check+1)=yaw(check+1)+2*pi;
    end
end
yaw=yaw-pi;


for i=1:nbDatab
Tb(:,:,i)=[rotz(yaw(i)) xb_repro(:,i); zeros(1,3) 1];
end

xdotb=computeDerivative(xb,dt);
xddotb=computeDerivative(xdotb,dt);



figure()
title('body trajectory')
axs=[-50 50 -50 50 -50 50];


for t=1:nbDatab
obs1=obs1_path(t,:)';
axis(axs)
hold on
plot3(xb_repro(1,:),xb_repro(2,:),xb_repro(3,:))
plot3(xb0(1),xb0(2),xb0(3),'ko')
plot3(xbg(1),xbg(2),xbg(3),'ro')
grid on
draw_sphere(obs1,po/2)
draw_sphere(obs2,po/2)

hold off
end




%% Link path defined and learned in joint space
global q_DMP qdot_DMP xl_repro xdotl_repro q_init q_fin
global wq Muq Sigmaq




q0=[0 -pi/2 0 0 0];
q_start=[1e-4 0 1e-4 pi/30 pi/30];
q_end=[1e-4 0 1e-4 pi/4 pi/3];



% still arms_flat is SE3(eye(4))
xl_start=arms_flat0{end}.fkine(qin(q0+q_start));
xl_end=arms_flat0{end}.fkine(qin(q0+q_end));
xl=mtraj(@lspb,xl_start.t',xl_end.t',nbDatal);
xl_nominal=xl;
xdotl=computeDerivative(xl',dt);
xddotl=computeDerivative(xdotl,dt);

[wl,Mul,Sigmal]=DMP_learning(nbDatal,3,alphal,dt, nbStatesl,Kl,Dl,xl');
[xl_repro,xdotl_repro]=DMP_repro(nbDatal,3,alphal,dt,nbStatesl,Kl,Dl ,wl,Mul,Sigmal,1,xl_start.t,xl_start.t,xl_end.t,0);

% link motion animation-did with not special reason

qs=zeros(nbDatal,dof);
for t=1:nbDatal     
    
    if t==1
        qs(t,:)=q0+q_start;
    else
        if t<=nbDatal % activating link 
            J=arms{dof}.end.jacob0(qin(prev_q));
            J=J(1:3,[1 5 7 9]);
            DMP_vel=(pinv(J,0.0001)*xdotl_repro(:,t));
            thetadot=DMP_vel;
            qs(t,:)=prev_q+[thetadot(1)';0;thetadot(2:end)]'*dt;
        else % no - activation link
            qs(t,:)=qs(t-1,:);
        end
    end
    prev_q=qs(t,:);
end


q_init=qs(1,:);
q_fin=qs(end,:);



[wq,Muq,Sigmaq]=DMP_learning(nbDatal,nbVarl,alphal,dt, nbStatesl,Kl,Dl,qs');

[q_DMP,qdot_DMP]=DMP_repro(nbDatal,nbVarl,alphal,dt,nbStatesl,Kl,Dl ,wq,Muq,Sigmaq,1,qs(1,:)',qs(1,:)',qs(end,:)',0);


figure()
title('link motion')
plot3(xl_repro(1,:),xl_repro(2,:),xl_repro(3,:),'k')
hold on
arms_flat0{end}.plot(qin(qs),'jointdiam',0.25,'jointcolor',[0.5 0.5 1],'tilesize',5,'nobase','noname','nowrist')


%% Data orgarnizing 
global xd yd zd 
xd=[ (1:200)' xb_repro(1,:)'];
yd=[(1:200)' xb_repro(2,:)' ];
zd=[ (1:200)' xb_repro(3,:)'];


%% slimulink here 
global iscol1 isrepro1 isaway1
iscol1=0; isrepro1=0; isaway1=0;







run('simulink_pd.mdl')

%% mapping 

rpy=interp1(output1(:,end),output1(:,4:6),1:nbDatab);
xyz=interp1(output1(:,end),output1(:,1:3),1:nbDatab);



% %% Flight Simulation 
% 
% % initialization 
% [xl_repro,xdotl_repro]=DMP_repro(nbDatal,3,alphal,dt,nbStatesl,Kl,Dl ,wl,Mul,Sigmal,1,xl_start.t,xl_start.t,xl_end.t,0);
% [q_DMP,qdot_DMP]=DMP_repro(nbDatal,nbVarl,alphal,dt,nbStatesl,Kl,Dl ,wq,Muq,Sigmaq,1,qs(1,:)',qs(1,:)',qs(end,:)',0);
% 
% q_sim=zeros(nbDatab,dof);
% 
% q_null1=zeros(nbDatab,dof);
% x0s=zeros(3,nbDatab);
% x0dots=zeros(3,nbDatab);
% figure()
% 
% 
% isaway1=0;
% iscol1=0;
% isrepro1=0;
% 
% % obs1=xb(:,nbDatab/4)+[0 1 -8]';
% 
% for t=1:nbDatab
%     
%      R=double(SE3(eye(3),xyz(t,:)))*rpy2tr(rpy(t,1),rpy(t,2),rpy(t,3));
%      obs1=obs1_path(t,:)';
%     
%     [arms,arms_flat]=Make_arm(LinkLength,R);
%     
%     if t==1
%         q_sim(t,:)=q0+q_start;
%     else
%             
%         J=arms{dof}.end.jacob0(qin(prev_q));
%         xe=arms{dof}.end.fkine(qin(prev_q));
%         xe=xe.t;
%         J=J(1:3,[1 5 7 9]);
%         [x01,J01,x0dot1,do1,nrst]=MK5(obs1,prev_q);
%         J01=J01(:,[1 3 4 5]);
%         
%             if do1<po
%                 iscol1=1;
%             end
%                  
%             if do1<po/3
%                  disp('colision-1')
%             end
% 
%                                    
%         qdot_null1=alphah(do1)*pinv(J01*(eye(dof-1)-pinv(J)*J),0.0001)*(alphao(do1)*x0dot1-J01*pinv(J)*xdotl_repro(:,t));
% %         disp(pinv(J01*(eye(dof-1)-pinv(J)*J),0.0001))
%            
%         
%         if ~iscol1
%         
%         qdot=qdot_DMP(:,t);    
% %         disp('following DMP- still didnt colide')
%         
%         
%         elseif iscol1 && ~isrepro1
%             if norm(pinv(J01*(eye(dof-1)-pinv(J)*J),0.0001))<1e-6
%                 qdot=pinv(J,0.001)*x0dot1;
%             else
%                 qdot=pinv(J,0.001)*xdotl_repro(:,t)+qdot_null1;
%             end
%                 
%             
%                 qdot=[qdot(1) ;0 ;qdot(2:end)];
% %         disp('following NULL')
%         
%         else 
%         qdot=qdot_DMP(:,t);    
% %         disp('following DMP- after colision')
%         
%         end
%         
%         isaway1=do1>po+0.5 && iscol1;
%                 
%         q_sim(t,:)=prev_q+qdot'*dt; 
%        
%        
%        
%        
%         if isaway1 && ~isrepro1
%            start_step=t;
%            start_pos=q_sim(t,:)';
%            [q_repro,qdot_repro]=DMP_repro(nbDatal,nbVarl,alphal,dt,nbStatesl,Kl,Dl ,wq,Muq,Sigmaq,start_step,start_pos,qs(1,:)',qs(end,:)',0);
%            qdot_DMP(:,t+1:end)=qdot_repro(:,t+1:end);
%            disp('----------------repro1-----------------------------')
%            isrepro1=1;
%         end
%         
%               
%     end
%     prev_q=q_sim(t,:);
% 
% 
%     arms_flat{end}.plot(qin(q_sim(t,:)),'nojoints','tilesize',10,'nobase','noname',...
%    'workspace',axs,'nowrist','fps',1000,'linkcolor','k')    
%     hold on
%     
% 
% 
% 
%     plot3(xb_repro(1,:),xb_repro(2,:),xb_repro(3,:))
%     plot3(xb0(1),xb0(2),xb0(3),'ko')
%     plot3(xbg(1),xbg(2),xbg(3),'ro')
%    
%     
% 
%     
%     if t~=1
%     quiver3(x01(1),x01(2),x01(3),x0dot1(1),x0dot1(2),x0dot1(3),'b')
%     end
% 
%     draw_sphere(obs1,po)
% end

%% data import 

rpy=interp1(output1(:,end)*10,output1(:,4:6),1:nbDatab);
xyz=interp1(output1(:,end)*10,output1(:,1:3),1:nbDatab);
q=interp1(output1(:,end)*10,output1(:,13:17),1:nbDatab);


%%
figure()

for t=1:nbDatab
    plot3(xb_repro(1,:),xb_repro(2,:),xb_repro(3,:))
    obs1=obs1_path(t,:)';
    R=double(SE3(eye(3),xyz(t,:)))*rpy2tr(rpy(t,1),rpy(t,2),rpy(t,3));
    hold on
    draw_drone(R,4,2,3)
       
%     x0=x0s(:,t);
%     x0dot=x0dots(:,t);
    axis([-50 50 -50 50 -25 25])

%     quiver3(x0(1),x0(2),x0(3),x0dot(1),x0dot(2),x0dot(3),'b')
%     
    draw_sphere(obs1,po/1.5)


    [arms,arms_flat]=Make_arm(LinkLength,R);
    

    arms_flat{end}.plot(qin([0 q(t,2:5)]),'nojoints','tilesize',10,'nobase','noname',...
   'workspace',axs,'nowrist','fps',100,'linkcolor','k','view',[-70 40],'workspace',[-50 50 -50 50 -25 25])    
    
    hold off

end



% 
% 
%             % trajectory reproduction 
%             
%             if do>po+1
%                 qdot_null=zeros(dof-1,1); 
%             else
%                start_step=t;
%                start_pos=arms{dof}.end.fkine(qin(prev_q)).t-xb(:,t);
%                [x_repro,xdot_repro]=DMP_repro(nbDatal,3,alphal,dt,nbStatesl,Kl,Dl,wl,Mul,Sigmal,start_step,start_pos,xl_start.t,xl_end.t);
%                xdotl(:,t+1:end)=xdot_repro(:,t+1:end);
%                isRepro=1;        
%             end
%        
        



















