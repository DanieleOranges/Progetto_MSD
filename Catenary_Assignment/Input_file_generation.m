%% Intro
clear
close all
clc

%% Points

A = [0;1.9];
B = [0;0];
C = [4;0.4];
D = [4;1.9];
E = [0.4;1.9];
F = [1.473;0.7];
G = [2;0.95];
H = [2.56;0.7];
I = [2.2;0.4];

%% Connections and Materials assignment
Connections = ["AD","BD","EG","FH","ID","IC"];
Materials   = ["Blue","Blue","Red","Red","Green","Green"];
n = 10*ones(length(Connections),1);

%% Element properties calculations
global EJ EA Area rho E_steel m

% Steel
E_steel   = 2.06*10^11;     % [N/m^2]
rho = 7800;                 % [Kg/m^3]   

t = [2,1.5,1.5]/1000;       % [m]
d = [60,40,25]/1000;        % [m]

J = pi/4.*[(d/2).^4-((d-2*t)/2).^4];
EJ          = E_steel*J;

Area        = pi.*((d/2).^2-((d-2*t)/2).^2);
EA          = E_steel*Area;
m           = rho*Area;

%% 

for k = 1 : length(Connections)
    
    P1 = eval(Connections{k}(1));
    P2 = eval(Connections{k}(2));
    [Lmax(k),nmin(k)] = element_size(P1,P2,Materials{k});
    [x(:,k),y(:,k)]       = mesh(P1,P2,n(k));
    
end

Ntot = length(x(:));
Node = 1:Ntot;
Elements = 1:(Ntot-length(Connections));

%% Nodes Assignment
delete('Fem_Assignment_input_prova')
fid = fopen('Fem_Assignment_input_prova','w');
fprintf(fid,'*NODES\n');
fprintf(fid,'%u    %f    %f\n',[Node;x(:)';y(:)']);
fprintf(fid,'*ENDNODES\n');

%% Elements Creation
fprintf(fid,'*BEAMS\n');

%%%% BISOGNA SELEZIONARE I NODI MASTER E SLAVE
% fprintf(fid,'%u    %u    %u    %u\n',[Elements;Node(2:end-1);Node(1:end)]);

fprintf(fid,'*ENDBEAMS\n');

%% Properties Creation
fprintf(fid,'*PROPERTIES\n');

fprintf(fid,'*ENDPROPERTIES\n');

% Close file
fclose(fid);
%% Functions

function [x,y] = mesh(P1,P2,n)

    x = [P1(1),P2(1)];
    y = [P1(2),P2(2)];
    
    x = interp1(1:length(x), x, linspace(1, length(x), n), 'linear');
    y = interp1(1:length(y), y, linspace(1, length(y), n), 'linear');
    
    figure(1); hold on; grid on;
    plot(x,y,'b-*')

end

function [Lmax,nmin] = element_size(P1,P2,Materials)

    global EJ m
    
    eta = 1.5; 
    OmegaMax = 20*(2*pi);      % [rad/s]
    
    switch Materials
        case 'Blue'
            index = 1;
        case 'Green'
            index = 2;
        case 'Red' 
            index = 3;
    end
    
    Lmax = sqrt(pi^2/eta/OmegaMax*sqrt(EJ(index)/m(index)));
    nmin = pdist([P1,P2])/Lmax;
    
end