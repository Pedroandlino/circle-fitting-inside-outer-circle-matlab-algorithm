clear all
clc
close all

% this is an example code developed by Pedro Lino %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%user inputs:
Rmax_user=8;
Rmin_user=4;

Max_iterations=1000;    %should be a high value, its the number of simulation cycles 
delta_t=0.01;           %should be a low value, in order to get small simulation steps
tolerance1=0.1;         %should be a low value, its related to a mean square error threshold
tolerance2=0.01;        %should be a low value, its related to a fine mean square error threshold
tolerance3=0.001;       %should be a low value, its related to a very small difference between 2 double values

%>>>plot flags<<<:
plot_flag=1;
final_plot_flag=0;


%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check user inputs:
if check_inputs(Rmax_user,Rmin_user,Max_iterations,delta_t,tolerance1,tolerance2,tolerance3)
    disp('Re-check your input values')
    return
end

%output:
Ncircles=0;

%constants:
Rmax=Rmax_user/Rmin_user;
Rmin=1;
x_pos=cos(linspace(-pi,pi));
y_pos=sin(linspace(-pi,pi));


%plot outer circle:
if plot_flag
    h_outer=setup_plot(Rmax_user,Rmin_user);
end

%plot boudary circle:
[H_all,pts,n_pts]=setup_field(Rmax,Rmin,Rmin_user,x_pos,y_pos,plot_flag);

%auxiliar cicle variables:
Ncircles_global=0;
pts_ini=pts;
remove_pts=[];
flag_done=0;

disp('Simulation Started, please stand by ...');

%decrement outlier points 
for losses=1:n_pts
    
    %simulation cycle
    for t=0:(Max_iterations) %note:could be while 1
        
        
%%%%%%%%%%%%%%%%%%%%%%% check for flag_done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %check if all circles are already inside
        Ncircles=0;
        k=0;
        for i=1:n_pts
            if ~ismember(i,remove_pts)
                k=k+1;
                if pts(i).pos*pts(i).pos'<=(Rmax-Rmin)^2
                    Ncircles=Ncircles+1;
                end
            end
        end
        if k==Ncircles
            %reached maximum already, becouse all of them fit inside!
            flag_done=1;
            break;
        end
        

        %special baundary/limit case!
        theta=2*pi/(k*2);
        Rmax_suggested=Rmin+Rmin/sin(theta);
        if Rmax_suggested-Rmax<tolerance3
            theta=(2*pi)/k;
            Ncircles=k;
            flag_done=1;
            
            k=0;
            for i=1:n_pts
                if ~ismember(i,remove_pts)
                    pts(i).pos(1)=(Rmax-Rmin)*cos(theta*k);
                    pts(i).pos(2)=(Rmax-Rmin)*sin(theta*k);

                    if plot_flag
                    %plot:
                        H_all(end-k).circ.XData=Rmin_user*( pts(i).pos(1) + x_pos);
                        H_all(end-k).circ.YData=Rmin_user*( pts(i).pos(2) + y_pos);
                    end
                    k=k+1;
                end
            end
            break;
        end
        
        %exact corners check at end of simulation cycle:
        if t==Max_iterations
            CM=[];
            the_error=[];
            for i=1:n_pts
                if ~ismember(i,remove_pts)
                    CM=[CM;pts(i).pos];
                end
            end
            CM=mean(CM,1);
            for i=1:n_pts
                if ~ismember(i,remove_pts)
                    pts(i).pos=pts(i).pos-CM;
                    pts(i).vel=[0 0];
                    
                    %if still inside, measure error
                    if pts(i).pos*pts(i).pos'>(Rmax-Rmin)^2
                        the_error=[the_error; abs(pts(i).pos*pts(i).pos' - (Rmax-Rmin)^2)];
                    end
                    if plot_flag
                        %plot
                        H_all(i).circ.XData = Rmin_user*(pts(i).pos(1) + x_pos);
                        H_all(i).circ.YData = Rmin_user*(pts(i).pos(2) + y_pos);
                    end
                end
            end
            
            if mean(the_error)<tolerance1
                %second minimization cycle for fine tunning:
                [the_fine_error,sugested_Ncircles]=fine_tunning(Rmin,Rmin_user,Rmax,H_all,Max_iterations,delta_t/100,pts,n_pts,remove_pts,plot_flag,x_pos,y_pos);
                
                if the_fine_error<tolerance2
                    disp(['Circles were considered to  be inside, with the mean square error: ' num2str(the_fine_error)])
                    Ncircles=sugested_Ncircles;
                    flag_done=1;
                    break;
                end
            end
        end
        
        
%%%%%%%%%%%%%%%%%%%%%%%% actual simulation cycle %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %transfer function
        if t<=Max_iterations
            Rnow=(Rmax+Rmin)*(1 -t/Max_iterations) + (Rmax)*(t/Max_iterations);
        end
        

        for i=1:n_pts
            if ~ismember(i,remove_pts)

                %squished by the red circle:
                if pts(i).pos*pts(i).pos'<(Rnow)^2 && (Rnow-Rmin)^2<pts(i).pos*pts(i).pos'
                    %converge:
                    pts(i).vel = pts(i).vel -pts(i).pos/norm(pts(i).pos);
                end

                %expected position after uptade:
                expected = pts(i).pos + delta_t*pts(i).vel;

                %inter collisions + best position to be now:
                for j=1:n_pts
                    if j~=i && ~ismember(j,remove_pts)

                        d = expected - pts(j).pos;

                        if d*d'<=(2*Rmin)^2
                            vel_proj=(pts(i).vel*d')*d/(d*d');
                            %newton collision:
                            pts(i).vel = pts(i).vel - vel_proj/2;
                            pts(j).vel = pts(j).vel + vel_proj/2;

                            %corrected pos:
                            expected = pts(j).pos + 2*Rmin*d/norm(d);

                        end
                    end
                end

                %settle
                pts(i).pos = expected;

                %clear outliers
                if pts(i).pos*pts(i).pos'>(Rnow)^2
                    pts(i).pos=[Inf,Inf];
                    pts(i).vel=[0,0];
                    remove_pts=[remove_pts, i];
                end 


                %plot:
                if plot_flag
                    H_all(i).circ.XData=Rmin_user*(pts(i).pos(1)+x_pos);
                    H_all(i).circ.YData=Rmin_user*(pts(i).pos(2)+y_pos);
                end
            end


        end

        %count
        Ncircles=0;
        for i=1:n_pts
            if ~ismember(i,remove_pts)
                 if pts(i).pos*pts(i).pos'<(Rmax-Rmin)^2
                     Ncircles=Ncircles+1;
                 end
            end
        end

        if plot_flag
            title(['Npts: ' num2str(n_pts) ' | i: ' num2str(losses+1) ' | N: ' num2str(Ncircles) ' | MAX: ' num2str(Ncircles_global) ]);

            %plot outer
            h_outer.XData=Rnow*Rmin_user*x_pos;
            h_outer.YData=Rnow*Rmin_user*y_pos;
            xlabel(['t: ' num2str(t)]);
            drawnow;
        end

    end
    
    %update max:
    if Ncircles_global<Ncircles
        Ncircles_global = Ncircles;
        if plot_flag
            title(['Npts: ' num2str(n_pts) ' | i: ' num2str(losses-1) ' | N= ' num2str(Ncircles) ' | MAX: ' num2str(Ncircles_global) ]);
        end
    end
    
    if flag_done
        if plot_flag
            ylabel('FLAG DONE');
        end
        break;
    end
    
    %reset and remove +1
    pts=pts_ini;
    remove_pts=[];
    for i=1:losses
        pts(i).pos=[Inf,Inf];
        pts(i).vel=[0 0];
        remove_pts=[remove_pts, i];
        if plot_flag
            H_all(i).circ.XData=Rmin_user*(pts(i).pos(1)+x_pos);
            H_all(i).circ.YData=Rmin_user*(pts(i).pos(2)+y_pos);
        end
    end
    
end


disp(['Ended! Max circ: ' num2str(Ncircles)]);

%final plot:
if final_plot_flag || plot_flag
    plot_final_result(pts,remove_pts,n_pts,Rmax_user,Rmin_user,Ncircles)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h=setup_plot(R,offset)
    h_p=figure(1);
    %h_p.Position=[1295 376 560 420];
    x_pos=R*cos(linspace(-pi,pi));
    y_pos=R*sin(linspace(-pi,pi));
    plot(x_pos,y_pos,'k--');hold on
    axis square
    xlim([-(R+offset) (R+offset)]);
    ylim([-(R+offset) (R+offset)]);
    x_pos=(R+offset)*cos(linspace(-pi,pi));
    y_pos=(R+offset)*sin(linspace(-pi,pi));
    h=plot(x_pos,y_pos,'r');
end

function [H,pts,n_pts]=setup_field(Rmax,Rmin,Rmin_user,x_pos,y_pos,plot_flag)
    H=struct;
    aux_pts=struct;
    pts=struct;
    cols=ceil(Rmax/double(Rmin));
    
    %fill with hexagonal space
    k=1;
    d=Rmin*realsqrt(3);
    N=ceil(Rmax/d);
    for i=-N:N
        
        if rem(i,2)==0
            for j=1:cols+1
                aux_pts(k).pos=[i*d,(j-1)*2*Rmin-(cols)*Rmin];
                %H(k).circ=plot(pts(k).pos(1)+x_pos,pts(k).pos(2)+y_pos,'b');
                k=k+1;
            end
        else
            for j=1:cols
                aux_pts(k).pos=[i*d,(j-1)*2*Rmin-(cols)*Rmin+Rmin];
                %H(k).circ=plot(pts(k).pos(1)+x_pos,pts(k).pos(2)+y_pos,'b');
                k=k+1;
            end
        end
    end
    
    n_aux_pts=k-1;
    k=1;
    %remove outliers
    for i=1:n_aux_pts
        if aux_pts(i).pos*aux_pts(i).pos'<(Rmax+Rmin)^2
            pts(k).pos=aux_pts(i).pos;
            pts(k).vel=[0 0];
            if plot_flag
                H(k).circ=plot(Rmin_user*(pts(k).pos(1)+x_pos),Rmin_user*(pts(k).pos(2)+y_pos),'b');
            end
            k=k+1;
        end
    end
    n_pts=k-1;
end


function plot_final_result(pts,remove_pts,n_pts,Rmax_user,Rmin_user,Ncircles)
    figure(1);
    hold off
    x_pos=Rmax_user*cos(linspace(-pi,pi));
    y_pos=Rmax_user*sin(linspace(-pi,pi));
    plot(x_pos,y_pos,'r');hold on
    axis square
    xlim([-Rmax_user Rmax_user]);
    ylim([-Rmax_user Rmax_user]);
    
    x_pos=Rmin_user*cos(linspace(-pi,pi));
    y_pos=Rmin_user*sin(linspace(-pi,pi));  
    
    for i=1:n_pts
        if ~ismember(i,remove_pts)
            plot(Rmin_user*pts(i).pos(1)+x_pos,Rmin_user*pts(i).pos(2)+y_pos,'b');
        end
    end
    
    title(['Rmax: ' num2str(Rmax_user) ' Rmin: ' num2str(Rmin_user) ' Max. Inside: ' num2str(Ncircles)]);
end

function [mean_error,sugested_Ncircles]=fine_tunning(Rmin,Rmin_user,Rnow,H_all,Max_iterations,delta_t,pts,n_pts,remove_pts,plot_flag,x_pos,y_pos)
    if plot_flag
        ylabel('fine tunning');
        color=H_all(1).circ.Color;
    end
    for t=0:Max_iterations
        for i=1:n_pts
            if ~ismember(i,remove_pts)

                %attactaion to center:
                pts(i).vel = pts(i).vel +([0 0]-pts(i).pos);

                expected = pts(i).pos + delta_t*pts(i).vel;

                %inter collisions, best position to be now:
                for j=1:n_pts
                    if j~=i && ~ismember(j,remove_pts)

                        d = expected - pts(j).pos;

                        if d*d'<=(2*Rmin)^2
                            vel_proj=(pts(i).vel*d')*d/(d*d');
                            %newton collision:
                            pts(i).vel = pts(i).vel - vel_proj/2;
                            pts(j).vel = pts(j).vel + vel_proj/2;

                            %corrected pos:
                            expected = pts(j).pos + 2*Rmin*d/norm(d);

                        end
                    end
                end


                pts(i).pos = expected;% + delta_t*pts(i).vel;

                %plot:
                if plot_flag
                    H_all(i).circ.XData=Rmin_user*(pts(i).pos(1)+x_pos);
                    H_all(i).circ.YData=Rmin_user*(pts(i).pos(2)+y_pos);
                    H_all(i).circ.Color=[0 1 1];
                end
            end
        end
        drawnow
    end
    
    %check new error
    the_error=[];
    sugested_Ncircles=0;
    for i=1:n_pts
        if ~ismember(i,remove_pts)
            %if still inside, measure error
            if pts(i).pos*pts(i).pos'>(Rnow-Rmin)^2
                the_error=[the_error; abs(pts(i).pos*pts(i).pos' - (Rnow-Rmin)^2)];
            end
            sugested_Ncircles=sugested_Ncircles+1;
            if plot_flag
                H_all(i).circ.Color=color;
            end
        end
    end
    
    if plot_flag
        ylabel('');
    end
    
    mean_error=mean(the_error);
end

function quit_by_error=check_inputs(Rmax_user,Rmin_user,Max_iterations,delta_t,tolerance1,tolerance2,tolerance3)
    quit_by_error=0;
    
    if Rmin_user<0 || Rmax_user<0
        quit_by_error=1;
        warning('Rmax and Rmin must be positive!');
    end
    
    if Rmax_user<Rmin_user
        quit_by_error=1;
        warning('Rmax inferior to Rmin!');
    end
    
    if Max_iterations<200
        quit_by_error=1;
        warning('Max_iterations should be more, to guaranty a good minimization cycle')
    end
    
    if delta_t>0.1
        quit_by_error=1;
        warning('delta_t should be smaller, to guaranty a good minimization cycle')
    end
    
    if tolerance1>0.1
        quit_by_error=1;
        warning('tolerance1 should be smaller, to guaranty precision')
    end
    
    if tolerance2>tolerance1
        quit_by_error=1;
        warning('tolerance2 should be smaller than tolerance1')
    end
    
    if tolerance3>0.01
        quit_by_error=1;
        warning('tolerance3 should be smaller, to guaranty precision comparison between doubles')
    end
end
