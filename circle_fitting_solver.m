clear all
clc
close all

% this is an example code developed by Pedro Lino %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%inputs:
Rmax=4;
Rmin=2;

Max_iterations=1000; %should be a high value, it won't necessarily reach it 
Max_crowns=10; %should be a high value, it won't necessarily reach it 
delta_t=0.5; %should be a low value, in order to get small steps

%>>>plot flag<<<:
plot_flag=0;


%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%protection:
if Rmax<2*Rmin || Max_iterations<100 || Max_crowns <1 || delta_t>0.5
    disp('Naughty boy!');
    return;
end


%output:
Ncircles=0;

%constants:


if plot_flag
    x_pos=Rmin*cos(linspace(-pi,pi));
    y_pos=Rmin*sin(linspace(-pi,pi));
else
    x_pos=[];
    y_pos=[];
end


%plot outer circle:
if plot_flag
    setup_plot(Rmax);
end

%plot inner crown:
[H_all,pts_in,n_pts_in]=setup_crown(Rmax,Rmin,x_pos,y_pos);

for crown=1:Max_crowns %note:could be while 1
    pts_out=pts_in;
    n_pts_out=n_pts_in;
    Ncircles=Ncircles+n_pts_out;
    [H,pts_in,n_pts_in]=gen_crown(Rmax-crown*(2*Rmin),Rmin,x_pos,y_pos,pts_out,n_pts_out,plot_flag);
    
    
    if n_pts_in<1
        %place one in center and see if fits!
        flag=0;
        for i=1:n_pts_out
            d=pts_out(i).pos;
            if d*d'<=(2*Rmin)^2
                %violation:
                flag=1;
                break;
            end
        end
        if flag==0
            %it fits, lets finish:
            if plot_flag
                H(1).circ=plot(x_pos,y_pos,'b');
                H_all=[H_all,H];
                Ncircles=Ncircles+1;
            end
        end
        break;
        
    end
    
    
    
    %simulation time:
    for t=0:Max_iterations %note:could be while 1
        count_ok=0;
        
        %physics collisions:
        for i=1:n_pts_in
            CM=[];
            
            %centrifuge:
            vel=pts_in(i).pos/norm(pts_in(i).pos);
            expected = pts_in(i).pos + delta_t*vel;
            
            
            %outer crown collisions:
            for j=1:n_pts_out
                expected = pts_in(i).pos + delta_t*vel;
                d=expected - pts_out(j).pos;
                
                if d*d'<=(2*Rmin)^2
                    vel_proj=(vel*d')*d/(d*d');
                    CM = pts_out(j).pos + 2*Rmin*d/norm(d);
                    pts_in(i).pos=CM;
                    vel = vel - vel_proj;
                end
                
            end
            
            %inter collisions:
            for j=1:n_pts_in
                if j~=i
                    expected = pts_in(i).pos + delta_t*vel;
                    d=expected - pts_in(j).pos;

                    if d*d'<=(2*Rmin)^2
                        vel_proj=(vel*d')*d/(d*d');
                        CM = pts_in(j).pos + 2*Rmin*d/norm(d);
                        pts_in(i).pos=CM;
                        vel = vel - vel_proj;
                    end
                    
                end
            end

            if isempty(CM)
                pts_in(i).pos=pts_in(i).pos + delta_t*vel;
            else
                %point stopped moving:
                count_ok=count_ok+1;
                pts_in(i).pos=CM;
            end
            
            
            %plot:
            if plot_flag
                H(i).circ.XData=pts_in(i).pos(1)+x_pos;
                H(i).circ.YData=pts_in(i).pos(2)+y_pos;
                xlabel([num2str(count_ok) '/' num2str(n_pts_in) ' - ' num2str(t)])
                drawnow;
            end
            
        end
        
        %break condition:
        if count_ok==n_pts_in
            break;
        end
    end
    
    
    
    %fix plot:
    if plot_flag
        title([num2str(crown+1) ' Crowns completed']);
        H_all=[H_all,H];
    end
    
end

if plot_flag
    title(['Ended! Max circ: ' num2str(Ncircles)]);
end
disp(['Ended! Max circ: ' num2str(Ncircles)]);

%nice exit:
clear d flag H i n_pts_in n_pts_out plot_flag pts_in pts_out ring x_pos y_pos 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setup_plot(R)
    h_p=figure(1);
    h_p.Position=[1295 376 560 420];
    x_pos=R*cos(linspace(-pi,pi));
    y_pos=R*sin(linspace(-pi,pi));
    plot(x_pos,y_pos,'r');hold on
    axis square
    xlim([-R R]);
    ylim([-R R]);
end

function [H,pts,n_pts]=setup_crown(Rmax,Rmin,x_pos,y_pos)
    H=struct;
    pts=struct;
    
    n_pts=floor(2*pi*(Rmax-Rmin)/(2*Rmin));
    
    theta=2*pi/double(n_pts);
    
    flag=0;
    for i=1:n_pts
        pts(i).pos=(Rmax-Rmin)*[cos((i-1)*theta),sin((i-1)*theta)];
        
        
        %check inter collision:
        for j=1:max(size(pts))
            if j~=i
                d=pts(i).pos-pts(j).pos;
                if d*d'<(2*Rmin)^2
                    flag=1;
                    break;
                end
            end
        end
        
        if ~isempty(x_pos) && flag==0
            H(i).circ=plot(pts(i).pos(1)+x_pos,pts(i).pos(2)+y_pos,'b');
        end
    end
    
    if flag
        %the the only other possible solution
        n_pts=1;
        pts(2:end)=[];
        H(2:end)=[];
        test_pos = (Rmax-Rmin)*[-1,0];
        d=test_pos-pts(1).pos;
        if d*d'>=(2*Rmin)^2
            pts(2).pos=test_pos;
            n_pts=n_pts+1;
            if ~isempty(x_pos)
                H(i).circ=plot(test_pos(1)+x_pos,test_pos(2)+y_pos,'b');
            end
        end
    end

end

function [H,pts,n_pts]=gen_crown(Rmax,Rmin,x_pos,y_pos,pts_in,n_pts_in,plot_flag)
    H=struct;
    pts=struct;
    
    %interlock:
    if n_pts_in>1
        cm=(pts_in(1).pos+pts_in(2).pos)/2;
        theta=atan2(cm(2),cm(1));
        Rot=[cos(theta) -sin(theta);sin(theta) cos(theta)];
    else
        Rot=eye(2);
    end
    
    n_pts=floor(2*pi*(Rmax-Rmin)/(2*Rmin));
    if n_pts<1
        return
    end
    theta=2*pi/double(n_pts);
    
        
    for i=1:n_pts
        pts(i).pos=(Rot*(Rmax-Rmin)*[cos((i-1)*theta);sin((i-1)*theta)])';

        if plot_flag
            H(i).circ=plot(pts(i).pos(1)+x_pos,pts(i).pos(2)+y_pos,'b');
        end
    end

end