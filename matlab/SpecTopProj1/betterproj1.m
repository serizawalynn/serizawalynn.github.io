%%Define ball
%mass of the ball
mass=1;

%radius
radius=8.5;

%gravity
gravity=-9.8;

%Initial Position of the center of the ball
ball_pos_xi=25;
ball_pos_yi=25;
ball_pos_zi=10;
ball_pos_center=[ball_pos_xi,ball_pos_yi,ball_pos_zi];

%Initial velocity of the ball
ball_vel_xi=0;
ball_vel_yi=0;
ball_vel_zi=0;
ball_vel_center=[ball_vel_xi,ball_vel_yi,ball_vel_zi];

%Initial acceleration of the ball
ball_acc_xi=0;
ball_acc_yi=0;
ball_acc_zi=gravity+0;
ball_acc_center=[ball_acc_xi,ball_acc_yi,ball_acc_zi];

%%define the net
%Initial length of each string in the net
initial_length=.95;

%Dimensions of the net
rows=51;
columns=51;

%spring stiffness
k=999;

%Build arrays for the positions, velocities, and acceleration of the nodes
net_pos_x=zeros(columns,rows);
net_pos_y=zeros(columns,rows);
net_pos_z=zeros(columns,rows);
net_vel_x=zeros(columns,rows);
net_vel_y=zeros(columns,rows);
net_vel_z=zeros(columns,rows);
net_acc_x=zeros(columns,rows);
net_acc_y=zeros(columns,rows);
net_acc_z=zeros(columns,rows);

%Iterate through the rows to dictate their initial positions
%We assume that the net starts at z=0
%If someone wants to place the net at some other z, the initial height of
%the ball can be adjusted instead.
counter_x=0;
counter_y=0;
for i=1:rows
    for j=1:columns
        net_pos_x(i,j)=(counter_x);
        net_pos_y(i,j)=(counter_y);
        net_pos_z(i,j)=0;
        counter_x=counter_x+1;
        if counter_x==columns
            counter_x=0;
            counter_y=counter_y+1;
        else
            continue
        end
    end
end

%%Code for movie/visuals
my_figure=figure(1);

axis tight
axis equal
v=VideoWriter("filename.mp4");
%video fps
v.FrameRate=500;
open(v)

%%Here, one can alter these variables to alter output of the code

%Timestep
timestep=.0001;

%Duration of the simulation
seconds=5;

%The angle of viewing in the xy-plane
azimuth=0;

%Angle of viewing between the z-axis and the xy-plane
elevation=0;

for time=0:timestep:seconds
    %%Plot sphere, net at current dt
    time
    dt=timestep;
    %Clear plot
    cla()
    
    %Turn on the gridlines
    axis equal
    
    %Dictate angle of viewing
    view(azimuth,elevation)
    
    %Plot the ball
    [X,Y,Z]=sphere;
    X=(X*radius)+ball_pos_center(1);
    Y=(Y*radius)+ball_pos_center(2);
    Z=(Z*radius)+ball_pos_center(3);
    plot3(X,Y,Z)
    %plot3(ball_pos_center(1),ball_pos_center(2),ball_pos_center(3),'-o','MarkerSize',15)
    
    %Keep the current plot
    hold on
   
    %Set the limits of the axes
    axis([-1 52 -1 52 -2 30])
    %surf(net_pos_x,net_pos_y,net_pos_z,'FaceAlpha',0.25)
    surf(net_pos_x,net_pos_y,net_pos_z)
    
    xlabel("X")
    ylabel("Y")
    zlabel("Z")
    
    %%%Calculate the next position of the ball
    %%Determine forces created due to the interaction between the ball and
    %%the net
    
    net_ball_vector=[0,0,0];
    oppo_vec=[0,0,0];
    for i=2:(rows-1)
        for j=2:(columns-1)
            net_pos_vec=[net_pos_x(i,j),net_pos_y(i,j),net_pos_z(i,j)];
            if (radius-norm(ball_pos_center-net_pos_vec))>0
                %Force from the ball
                L=norm(ball_pos_center-net_pos_vec);
                vector_mag=(k)*(radius-L);
                vector_direction=-(net_pos_vec-ball_pos_center)/L;
                vector=vector_mag*vector_direction;
                net_ball_vector=net_ball_vector+vector;
                oppo_vec=oppo_vec+(((k)*(radius-L)*(net_pos_vec-ball_pos_center))/L);
                
                %%Force from the other nodes
                net_acc_x(i,j)=oppo_vec(1,1)/mass;
                net_acc_y(i,j)=oppo_vec(1,2)/mass;
                net_acc_z(i,j)=oppo_vec(1,3)/mass;
                Tension_up=tension(net_pos_x(i-1,j),net_pos_y(i-1,j),net_pos_z(i-1,j),net_pos_x(i,j),net_pos_y(i,j),net_pos_z(i,j),k,initial_length);
                Tension_down=tension(net_pos_x(i+1,j),net_pos_y(i+1,j),net_pos_z(i+1,j),net_pos_x(i,j),net_pos_y(i,j),net_pos_z(i,j),k,initial_length);
                Tension_left=tension(net_pos_x(i,j-1),net_pos_y(i,j-1),net_pos_z(i,j-1),net_pos_x(i,j),net_pos_y(i,j),net_pos_z(i,j),k,initial_length);
                Tension_right=tension(net_pos_x(i,j+1),net_pos_y(i,j+1),net_pos_z(i,j+1),net_pos_x(i,j),net_pos_y(i,j),net_pos_z(i,j),k,initial_length);
                Tension=((Tension_up+Tension_down+Tension_left+Tension_right)/mass);
                
                %Euler Method on node
                net_acc_x(i,j)=net_acc_x(i,j)+Tension(1,1);
                net_acc_y(i,j)=net_acc_y(i,j)+Tension(1,2);
                net_acc_z(i,j)=net_acc_z(i,j)+Tension(1,3);
                net_vel_x(i,j)=net_vel_x(i,j)+(net_acc_x(i,j)*dt);
                net_vel_y(i,j)=net_vel_y(i,j)+(net_acc_y(i,j)*dt);
                net_vel_z(i,j)=net_vel_z(i,j)+(net_acc_z(i,j)*dt);
                net_pos_x(i,j)=net_pos_x(i,j)+(net_vel_x(i,j)*dt)+(0.5*(dt^2)*net_acc_x(i,j));
                net_pos_y(i,j)=net_pos_y(i,j)+(net_vel_y(i,j)*dt)+(0.5*(dt^2)*net_acc_y(i,j));
                net_pos_z(i,j)=net_pos_z(i,j)+(net_vel_z(i,j)*dt)+(0.5*(dt^2)*net_acc_z(i,j));
               
            else
                %Calculate net tension exerted on node
                Tension_up=tension(net_pos_x(i-1,j),net_pos_y(i-1,j),net_pos_z(i-1,j),net_pos_x(i,j),net_pos_y(i,j),net_pos_z(i,j),k,initial_length);
                Tension_down=tension(net_pos_x(i+1,j),net_pos_y(i+1,j),net_pos_z(i+1,j),net_pos_x(i,j),net_pos_y(i,j),net_pos_z(i,j),k,initial_length);
                Tension_left=tension(net_pos_x(i,j-1),net_pos_y(i,j-1),net_pos_z(i,j-1),net_pos_x(i,j),net_pos_y(i,j),net_pos_z(i,j),k,initial_length);
                Tension_right=tension(net_pos_x(i,j+1),net_pos_y(i,j+1),net_pos_z(i,j+1),net_pos_x(i,j),net_pos_y(i,j),net_pos_z(i,j),k,initial_length);
                Tension=((Tension_up+Tension_down+Tension_left+Tension_right)/mass);
                
                %Euler method on net
                net_acc_x(i,j)=Tension(1,1);
                net_acc_y(i,j)=Tension(1,2);
                net_acc_z(i,j)=Tension(1,3);
                net_vel_x(i,j)=net_vel_x(i,j)+(net_acc_x(i,j)*dt);
                net_vel_y(i,j)=net_vel_y(i,j)+(net_acc_y(i,j)*dt);
                net_vel_z(i,j)=net_vel_z(i,j)+(net_acc_z(i,j)*dt);
                net_pos_x(i,j)=net_pos_x(i,j)+(net_vel_x(i,j)*dt)+(0.5*(dt^2)*net_acc_x(i,j));
                net_pos_y(i,j)=net_pos_y(i,j)+(net_vel_y(i,j)*dt)+(0.5*(dt^2)*net_acc_y(i,j));
                net_pos_z(i,j)=net_pos_z(i,j)+(net_vel_z(i,j)*dt)+(0.5*(dt^2)*net_acc_z(i,j));
            end
        end
    end
   
   %%Calculate the new position,velocity of the ball
    ball_acc_center=((net_ball_vector)./mass)+[0 0 -9.8];
    
    ball_vel_center=ball_vel_center+(dt*ball_acc_center);

    ball_pos_center=ball_pos_center+(ball_vel_center*dt)+(.5*(dt^2)*ball_acc_center);
    
    %Code for visuals
    Mv=getframe(gcf);
    writeVideo(v,Mv)
    
end
    
close(v)

%Function to determine direction of the tension vector
function [direction_vec] = direction(x_pulling,y_pulling,x_node,y_node)
if (x_pulling-x_node>0) && (y_pulling-y_node)>0
    direction_vec=[-1 -1];   
elseif (x_pulling-x_node)>0 && (y_pulling-y_node)==0
    direction_vec=[-1 0];
elseif x_pulling-x_node>0 && y_pulling-y_node<0
    direction_vec=[-1 1];
elseif x_pulling-x_node==0 && y_pulling-y_node>0
    direction_vec=[0 -1];
elseif x_pulling-x_node==0 && y_pulling-y_node==0
    direction_vec=[0 0];
elseif x_pulling-x_node==0 && y_pulling-y_node<0
    direction_vec=[0 1];
elseif x_pulling-x_node<0 && y_pulling-y_node>0
    direction_vec=[1 -1];
elseif x_pulling-x_node<0 && y_pulling-y_node==0
    direction_vec=[1 0];
elseif x_pulling-x_node<0 && y_pulling-y_node<0
    direction_vec=[1 1];
end
end
%Function to calculate tension
function [tension_vector] = tension(x_pulling,y_pulling,z_pulling,np_x,np_y,np_z,k,initial_length)
pull=[x_pulling,y_pulling,z_pulling];
node=[np_x,np_y,np_z];
change_x=abs(pull(1,1)-node(1,1));
change_y=abs((pull(1,2)-node(1,2)));
change_z=((pull(1,3)-node(1,3)));
net_change_spring=initial_length-sqrt((change_x^2)+(change_y^2)+(change_z^2));
T=k*net_change_spring;
phi=atan(change_z/sqrt((change_x^2)+(change_y^2)));
Txy=T*cos(phi);
Tz=T*sin(phi);
theta=atan(change_y/change_x);
Tx=Txy*cos(theta);
Ty=Txy*sin(theta);
direction_vec=direction(pull(1,1),pull(1,2),node(1,1),node(1,2));
tension_vector=[direction_vec(1,1)*Tx,direction_vec(1,2)*Ty,-Tz];
end