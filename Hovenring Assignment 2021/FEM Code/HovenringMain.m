%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical System Dynamics
% FEM script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
close all; 
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load input file and assemble structure
[file_i,xy,nnod,sizee,idb,ndof,incid,l,gamma,m,EA,EJ,T,posit,nbeam,pr]=loadstructure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw structure
dis_stru(posit,l,gamma,xy,pr,idb,ndof);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble mass and stiffness matricies
[M,K] = assem(incid,l,m,EA,EJ,T,gamma,idb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add concentrated elements
% Torsional spring - deck side
kT1 = 9e6;
i_ndof_spring = idb(1,3);
K(i_ndof_spring,i_ndof_spring) = K(i_ndof_spring,i_ndof_spring)+kT1;

% Torisional spring - tower side
kT2 = 9e6;
i_ndof_spring = idb(230,3);
K(i_ndof_spring,i_ndof_spring) = K(i_ndof_spring,i_ndof_spring)+kT2;

% TMD 1 - 04R36 - Mass and stiffness
m_tmd1 = 4.93;
k_tmd1 = 2.12E+03;
i_ndof1 = idb(149,1);
i_ndof2 = idb(231,1);
i_dof_tmd1 = [i_ndof1 i_ndof2];

M(i_ndof2,i_ndof2) = M(i_ndof2,i_ndof2) + m_tmd1;

K_tmd1 = [k_tmd1 -k_tmd1; -k_tmd1 k_tmd1];
K(i_dof_tmd1,i_dof_tmd1) = K(i_dof_tmd1,i_dof_tmd1)+K_tmd1;

% TMD 2 - 04R36 - Mass and stiffness
m_tmd1 = 5.1;
k_tmd1 = 888;
i_ndof1 = idb(149,1);
i_ndof2 = idb(232,1);
i_dof_tmd1 = [i_ndof1 i_ndof2];

M(i_ndof2,i_ndof2) = M(i_ndof2,i_ndof2) + m_tmd1;

K_tmd1 = [k_tmd1 -k_tmd1; -k_tmd1 k_tmd1];
K(i_dof_tmd1,i_dof_tmd1) = K(i_dof_tmd1,i_dof_tmd1)+K_tmd1;


% TMD 1 - 4RZ11 - Mass and stiffness
m_tmd1 = 5;
k_tmd1 = 4.00E+03;
i_ndof1 = idb(164,1);
i_ndof2 = idb(233,1);
i_dof_tmd1 = [i_ndof1 i_ndof2];

M(i_ndof2,i_ndof2) = M(i_ndof2,i_ndof2) + m_tmd1;

K_tmd1 = [k_tmd1 -k_tmd1; -k_tmd1 k_tmd1];
K(i_dof_tmd1,i_dof_tmd1) = K(i_dof_tmd1,i_dof_tmd1)+K_tmd1;

% TMD 2 - 4RZ11 - Mass and stiffness
m_tmd1 = 7;
k_tmd1 = 1.17E+04;
i_ndof1 = idb(164,1);
i_ndof2 = idb(234,1);
i_dof_tmd1 = [i_ndof1 i_ndof2];

M(i_ndof2,i_ndof2) = M(i_ndof2,i_ndof2) + m_tmd1;

K_tmd1 = [k_tmd1 -k_tmd1; -k_tmd1 k_tmd1];
K(i_dof_tmd1,i_dof_tmd1) = K(i_dof_tmd1,i_dof_tmd1)+K_tmd1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute natural frequencies and mode shapes
MFF = M(1:ndof,1:ndof);
MCF = M(ndof+1:end,1:ndof);
MFC = M(1:ndof,ndof+1:end);
MCC = M(ndof+1:end,ndof+1:end);

KFF = K(1:ndof,1:ndof);
KCF = K(ndof+1:end,1:ndof);
KFC = K(1:ndof,ndof+1:end);
KCC = K(ndof+1:end,ndof+1:end);

[modes, omega] = eig(MFF\KFF);
omega = sqrt(diag(omega));
[omega,i_omega] = sort(omega);
freq0 = omega/2/pi;
modes = modes(:,i_omega);

nmodes = 2;
scale_factor = 30;
for ii=1:nmodes
    mode = modes(:,ii);
    figure();
    diseg2(mode,scale_factor,incid,l,gamma,posit,idb,xy)
    title(['mode ',num2str(ii),' freq ',num2str(freq0(ii)),' Hz']);
end