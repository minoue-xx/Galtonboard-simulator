% Galtonboard Simulator
% 22,Oct,2021 by MKT@KMKTo (Twitter)
%
% Assumtions
% - Particles are point mass.
% - No interactions in between particles.
% - Particles contact pegs at the outer radius.
% - Particles bounce off pegs to normal force direction with restitution (e).

clear
close all

%% Initial setup
N_step = 10000; %Maximum steps
dt     = 0.01;  %[sec]
g      = -9.8;  %m/ss, gravity
e      = 0.4;%  %Coefficient of restitution
N_part = 100;  %Number of particles

N_layers = 8; % Number of layers. Even number recommemed
N_pegs_per_row = 20; %Total number of pegs per row
N_obj  = (N_layers+1)*N_pegs_per_row;   %Total number of pegs

r_obj  = 0.25;  %[m] Diameter of pegs
rng(0);         %Random initialization

%% Creation of video object
writerObj = VideoWriter('simulation.avi');
open(writerObj);

%% Particle initialization
nlayer = N_obj/N_pegs_per_row;
h_drop = nlayer+2;
step = 1:N_step;
for part = 1:N_part
    p(part).v = zeros(N_step,2);             %Speed (x,y direction), allocation
    p(part).x = zeros(N_step,2);             %Position (x,y direction), allocation
    p(part).x(1,1) =  randn*0.01;    %Initial position (x direction)
    p(part).x(1,2) =  h_drop;            %Initial position (y direction)
    p(part).inobj     = 0;              %Identifier for in-out particle
end

%% Peg definition
N_draw = 10;
xcircle = r_obj*sin([0:N_draw]/N_draw*2*pi); %For drawing circle radius
ycircle = r_obj*cos([0:N_draw]/N_draw*2*pi); %For drawing circle radius

for nobj = 1:N_obj
    ob(nobj).x(2) =  floor(nobj/N_pegs_per_row);     %Y position of the center of the peg
    ob(nobj).x(1) =  1*((mod(nobj,N_pegs_per_row)-N_pegs_per_row/2)+mod(ob(nobj).x(2),2)*0.5);     %X position of the center of the peg
    ox(nobj).xplt = ob(nobj).x(1)+xcircle; %For drawing circle radius
    ox(nobj).yplt = ob(nobj).x(2)+ycircle; %For drawing circle radius
end

%% Drawing of the pegs
screanW = 950;
screanH = 600;

figure('Position',[0 0 screanW screanH])
for nobj = 1:N_obj
    hold on
    plot(ox(nobj).xplt,ox(nobj).yplt,'Color',[0 0.4470 0.7410])
end
text(-N_pegs_per_row/2+2,(N_pegs_per_row-2)*screanH/screanW-2,['e=' num2str(e)]); %Coefficient of restitution

%% Start of simulation steps
view = 0;
for step = 2:N_step %Simulation steps

    if mod(step,5) == 0 && step < 600
        for part = length(p)+1:length(p)+10
            p(part).v = zeros(N_step,2);             %Speed (x,y direction), allocation
            p(part).x = zeros(N_step,2);             %Position (x,y direction), allocation
            p(part).x(1:step,1) =  randn*0.01;    %Initial position (x direction)
            p(part).x(1:step,2) =  h_drop;            %Initial position (y direction)
            p(part).inobj     = 0;              %Identifier for in-out particle
        end
        N_part = N_part + 10;
    end

    for part = 1:N_part % Time update for each particles

        if p(part).x(step-1,2)<=0  % If the particle reaches the floor, stop the integration
            p(part).x(step,2)=0;
            p(part).x(step,1)=p(part).x(step-1,1); %Stop movement of the particle
            x(part)=p(part).x(step,1);             %Record for drawding
            y(part)=p(part).x(step,2);             %Record for drawding
        else  % If not reaches the floor
            % Time integration with Euler method
            p(part).v(step,1) =  p(part).v(step-1,1);
            p(part).v(step,2) =  p(part).v(step-1,2)+g*dt;
            p(part).x(step,1) =  p(part).x(step-1,1)+p(part).v(step,1)*dt;
            p(part).x(step,2) =  p(part).x(step-1,2)+p(part).v(step,2)*dt;
            x(part)=p(part).x(step,1);             %Record for drawding
            y(part)=p(part).x(step,2);             %Record for drawding

            % Reflection of the particles toward pegs
            for nobj = 1:N_obj
                u  = p(part).x(step,:)-ob(nobj).x; %vector between particle vs peg
                nu = norm(u); %distance
                if nu<r_obj   %if distance is close enough
                    if p(part).inobj ~= 1   %if the particle is out of the peg area in the previous step
                        u = u/nu; %Normalized normal force direction
                        ut = (u*p(part).v(step,:)')*u; %Velocity component toward normal force direction
                        un = p(part).v(step,:)-ut; %Velocity component orthogonal to normal force directio
                        vupate = un-e*ut; % Update of the speed; flip the vector only to normal force direction
                        p(part).v(step,:)=vupate;
                        p(part).x(step,:)=ob(nobj).x+u*r_obj; %Bring the particle surface of the peg
                        p(part).inobj     = 1;
                    end
                else
                    p(part).inobj     = 0;
                end
            end
        end
    end

    if view == 1 % Number of skips for the visualization, 1: visualize every steps

        if exist('p1') % If already particles are drawn in the previous steps
            delete(p1) % Delete the previous particle
        end

        p1 = plot(x,y,'.','MarkerSize',6,'Color',[0.8500 0.3250 0.0980]); % Draw particles


        xlim([-N_pegs_per_row/2+1,N_pegs_per_row/2-1])
        ylim([0,(N_pegs_per_row-2)*screanH/screanW])
        drawnow

        % Recording to video
            frame = getframe(gcf);
            writeVideo(writerObj, frame);
        view = 0;
    end
    view = view + 1;

    xfloor = x(y<=0);
    bin=[-14:14];
    [N,edges] = histcounts(xfloor,bin);
    xpos = edges(1:end-1)+0.5;
    if ~exist('hbar','var')
        hbar = bar(xpos,N/N_part*(N_pegs_per_row-2)*screanH/screanW*3);
        hbar.FaceAlpha = 0.6;
    else
        hbar.XData = xpos;
        hbar.YData = N/N_part*(N_pegs_per_row-2)*screanH/screanW*3;
    end


    % If all the particles are on the floor, then stops the steps
    if max(y)<=0
        break
    end
end
% End of the loop

%% Visualization of histogram
bin=[-14:14];
[N,edges] = histcounts(x,bin);
xpos = edges(1:end-1)+0.5;
hold on
hbar = bar(xpos,N/N_part*(N_pegs_per_row-2)*screanH/screanW*3);   %Histogram
hbar.FaceAlpha = 0.6;
text(-N_pegs_per_row/2+2,(N_pegs_per_row-2)*screanH/screanW-3,['std=' num2str(std(x))]);
xnorm = [-14:0.1:14];
ynorm = normpdf(xnorm,0,std(x))*(N_pegs_per_row-2)*screanH/screanW*3;
plot(xnorm,ynorm,'-r','LineWidth',1.5);     %Normal distribution

% Recording the histogram to the video
frame = getframe(gcf);
for roop = 1:200
writeVideo(writerObj, frame);
end
close(writerObj)

%% Drawing the particle trajectory
figure('Position',[0 0 950 600])

% drawing pegs
for nobj = 1:N_obj
    hold on
    plot(ox(nobj).xplt,ox(nobj).yplt,'Color',[0 0.4470 0.7410])
end
% drawing trajectory
for part = 1:min(50,N_part)%N_part
    plot(p(part).x(1:step,1),p(part).x(1:step,2),'-b');
end

xlim([-N_pegs_per_row/2+1,N_pegs_per_row/2-1])
ylim([0,(N_pegs_per_row-2)*screanH/screanW])

