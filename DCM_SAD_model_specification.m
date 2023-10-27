   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DYNAMIC CAUSAL MODELLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear DCM

basepath = '/data/netapp01/work/rebekah/Paper_2/Data/1_HC/';
basepath2 = '/data/netapp01/work/rebekah/Paper_2/Data/2_SAD/';
basepath1 = '/data/netapp01/work/rebekah/Paper_2/Data/3_SAD_MDD/';


id = {'H005'
'H007'
'H012'
'H014'
'H022'
'H025'
'H030'
'H031'
'H033'
'H038'
'H042'
'H043'
'H044'
'H046'
'H047'
'H056'
'H058'
'H063'
'H064'
'H066'
'H067'
'H068'
'H071'
'H073'
'H078'
'H085'
'H093'
'H094'
'H095'
'H097'
'H098'
'H100'
'H101'
'H104'
'H108'
'H109'
'H116'
'H119'
'H121'
'H123'
'H125'
'H126'
'H131'
'H136'
'H138'
'H139'
'H140'
'H141'
'H149'
'H150'
'H153'
'H155'
'H157'
'H160'
'H162'
'H164'
'H165'
'H168'
'H170'
'H171'
'H172'
'H175'
'H176'
'H180'
'H187'
'H188'
'H193'
'H197'
'H198'
'H201'
'H202'
'H212'
'C001'
'C002'
'C004'
'C005'
'C006'
'C008'
'C009'
'C010'
'C011'
'C012'
'C015'
'C016'
'C017'
'C018'
'C019'
'C020'
'C022'
'C024'
'C101'
'C103'
'C104'
'C500'
'C502'
'C504'
'C505'
'S043'
'S049'
'S074'
'S080'
'S098'
'S101'
'S226'
'S228'
'S237'
'S246'
'S253'
'S501'
'S518'};

for i = 1:72
    curr_id = id(i);
 
    subject_num = i;
    
      
   ne_id = strcat(basepath, curr_id);  
   new_id = strcat(ne_id, '/1st_level/1st_level_additional'); 

   sav_id = strcat(new_id, '/DCM_Files_57_seperate_test');
   
   data_path = char(new_id);
   
   save_path = char(sav_id);
     
   cd(data_path)
   delete('DCM_Files_57_seperate_test\*')
   mkdir('DCM_Files_57_seperate_test')
  
% SPECIFICATION DCM "attentional modulation of backward connection"
%--------------------------------------------------------------------------
% To specify a DCM, you might want to create a template one using the GUI
% then use spm_dcm_U.m and spm_dcm_voi.m to insert new inputs and new
% regions. The following code creates a DCM file from scratch, which
% involves some technical subtleties and a deeper knowledge of the DCM
% structure.

load(fullfile(data_path,'SPM.mat'));

load(fullfile(data_path,'VOI_MPFC_1.mat'),'xY');
DCM.xY(1) = xY;
load(fullfile(data_path,'VOI_PCC_1.mat'),'xY');
DCM.xY(2) = xY;
load(fullfile(data_path,'VOI_left IPL_1.mat'),'xY');
DCM.xY(3) = xY;


DCM.n = length(DCM.xY);      % number of regions
DCM.v = length(DCM.xY(1).u); % number of time points

% Time series
%

DCM.Y.dt  = SPM.xY.RT;
DCM.Y.X0  = DCM.xY(1).X0;
for i = 1:DCM.n
    DCM.Y.y(:,i)  = DCM.xY(i).u;
    DCM.Y.name{i} = DCM.xY(i).name;
end

DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v);

% Experimental inputs

DCM.U.dt   =  SPM.Sess.U(1).dt;
DCM.U.name = [SPM.Sess.U.name];
DCM.U.u    = [SPM.Sess.U(1).u(33:end,1) ...
              SPM.Sess.U(2).u(33:end,1) ...
              SPM.Sess.U(3).u(33:end,1) ...
              SPM.Sess.U(4).u(33:end,1) ...
              SPM.Sess.U(5).u(33:end,1)];

DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
DCM.TE     = 0.035;

DCM.options.nonlinear  = 0;
DCM.options.two_state  = 0;
DCM.options.stochastic = 0;
DCM.options.centre    = 0;
DCM.options.induced    = 0;

% SPECIFICATION DCM 1 - "1x1"
%-------------------------------------------------------------------------
DCM.a = [1,1,1; 1,1,1; 1,1,1];
DCM.b = zeros(3,3,5); DCM.b(2,1,2)= 1; DCM.b(1,2,2)= 1; DCM.b(3,1,2)= 1; DCM.b(1,3,2)= 1; DCM.b(3,2,2)= 1; DCM.b(2,3,2)= 1; DCM.b(2,1,4)= 1; DCM.b(1,2,4)= 1; DCM.b(3,1,4)= 1; DCM.b(1,3,4)= 1; DCM.b(3,2,4)= 1; DCM.b(2,3,4)= 1;
DCM.c = [0 0 0 0 0; 1 1 0 1 0; 0 0 0 0 0];
DCM.d = zeros(3,3,0);

save(fullfile(save_path,'DCM_model_1_1x1.mat'),'DCM');






clear DCM
end


for i = 73:97
    curr_id = id(i);
 
    subject_num = i;
    
      
   ne_id = strcat(basepath1, curr_id);  
   new_id = strcat(ne_id, '/1st_level/1st_level_additional'); 

   sav_id = strcat(new_id, '/DCM_Files_57_seperate_test');
   
   data_path = char(new_id);
   
   save_path = char(sav_id);
     
   cd(data_path)
   delete('DCM_Files_57_seperate_test\*')
   mkdir('DCM_Files_57_seperate_test')
  
% SPECIFICATION DCM "attentional modulation of backward connection"
%--------------------------------------------------------------------------
% To specify a DCM, you might want to create a template one using the GUI
% then use spm_dcm_U.m and spm_dcm_voi.m to insert new inputs and new
% regions. The following code creates a DCM file from scratch, which
% involves some technical subtleties and a deeper knowledge of the DCM
% structure.

load(fullfile(data_path,'SPM.mat'));

load(fullfile(data_path,'VOI_MPFC_1.mat'),'xY');
DCM.xY(1) = xY;
load(fullfile(data_path,'VOI_PCC_1.mat'),'xY');
DCM.xY(2) = xY;
load(fullfile(data_path,'VOI_left IPL_1.mat'),'xY');
DCM.xY(3) = xY;


DCM.n = length(DCM.xY);      % number of regions
DCM.v = length(DCM.xY(1).u); % number of time points

% Time series
%

DCM.Y.dt  = SPM.xY.RT;
DCM.Y.X0  = DCM.xY(1).X0;
for i = 1:DCM.n
    DCM.Y.y(:,i)  = DCM.xY(i).u;
    DCM.Y.name{i} = DCM.xY(i).name;
end

DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v);

% Experimental inputs

DCM.U.dt   =  SPM.Sess.U(1).dt;
DCM.U.name = [SPM.Sess.U.name];
DCM.U.u    = [SPM.Sess.U(1).u(33:end,1) ...
              SPM.Sess.U(2).u(33:end,1) ...
              SPM.Sess.U(3).u(33:end,1) ...
              SPM.Sess.U(4).u(33:end,1) ...
              SPM.Sess.U(5).u(33:end,1)];

DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
DCM.TE     = 0.035;

DCM.options.nonlinear  = 0;
DCM.options.two_state  = 0;
DCM.options.stochastic = 0;
DCM.options.centre    = 0;
DCM.options.induced    = 0;

% SPECIFICATION DCM 1 - "1x1"
%-------------------------------------------------------------------------
DCM.a = [1,1,1; 1,1,1; 1,1,1];
DCM.b = zeros(3,3,5); DCM.b(2,1,2)= 1; DCM.b(1,2,2)= 1; DCM.b(3,1,2)= 1; DCM.b(1,3,2)= 1; DCM.b(3,2,2)= 1; DCM.b(2,3,2)= 1; DCM.b(2,1,4)= 1; DCM.b(1,2,4)= 1; DCM.b(3,1,4)= 1; DCM.b(1,3,4)= 1; DCM.b(3,2,4)= 1; DCM.b(2,3,4)= 1;
DCM.c = [0 0 0 0 0; 1 1 0 1 0; 0 0 0 0 0];
DCM.d = zeros(3,3,0);

save(fullfile(save_path,'DCM_model_1_1x1.mat'),'DCM');






clear DCM
end

for i = 98:110
    curr_id = id(i);
 
    subject_num = i;
    
      
   ne_id = strcat(basepath2, curr_id);  
   new_id = strcat(ne_id, '/1st_level/1st_level_additional'); 

   sav_id = strcat(new_id, '/DCM_Files_57_seperate_test');
   
   data_path = char(new_id);
   
   save_path = char(sav_id);
     
   cd(data_path)
   delete('DCM_Files_57_seperate_test\*')
   mkdir('DCM_Files_57_seperate_test')
  
% SPECIFICATION DCM "attentional modulation of backward connection"
%--------------------------------------------------------------------------
% To specify a DCM, you might want to create a template one using the GUI
% then use spm_dcm_U.m and spm_dcm_voi.m to insert new inputs and new
% regions. The following code creates a DCM file from scratch, which
% involves some technical subtleties and a deeper knowledge of the DCM
% structure.

load(fullfile(data_path,'SPM.mat'));

load(fullfile(data_path,'VOI_MPFC_1.mat'),'xY');
DCM.xY(1) = xY;
load(fullfile(data_path,'VOI_PCC_1.mat'),'xY');
DCM.xY(2) = xY;
load(fullfile(data_path,'VOI_left IPL_1.mat'),'xY');
DCM.xY(3) = xY;


DCM.n = length(DCM.xY);      % number of regions
DCM.v = length(DCM.xY(1).u); % number of time points

% Time series
%

DCM.Y.dt  = SPM.xY.RT;
DCM.Y.X0  = DCM.xY(1).X0;
for i = 1:DCM.n
    DCM.Y.y(:,i)  = DCM.xY(i).u;
    DCM.Y.name{i} = DCM.xY(i).name;
end

DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v);

% Experimental inputs

DCM.U.dt   =  SPM.Sess.U(1).dt;
DCM.U.name = [SPM.Sess.U.name];
DCM.U.u    = [SPM.Sess.U(1).u(33:end,1) ...
              SPM.Sess.U(2).u(33:end,1) ...
              SPM.Sess.U(3).u(33:end,1) ...
              SPM.Sess.U(4).u(33:end,1) ...
              SPM.Sess.U(5).u(33:end,1)];

DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
DCM.TE     = 0.035;

DCM.options.nonlinear  = 0;
DCM.options.two_state  = 0;
DCM.options.stochastic = 0;
DCM.options.centre    = 0;
DCM.options.induced    = 0;

% SPECIFICATION DCM 1 - "1x1"
%-------------------------------------------------------------------------
DCM.a = [1,1,1; 1,1,1; 1,1,1];
DCM.b = zeros(3,3,5); DCM.b(2,1,2)= 1; DCM.b(1,2,2)= 1; DCM.b(3,1,2)= 1; DCM.b(1,3,2)= 1; DCM.b(3,2,2)= 1; DCM.b(2,3,2)= 1; DCM.b(2,1,4)= 1; DCM.b(1,2,4)= 1; DCM.b(3,1,4)= 1; DCM.b(1,3,4)= 1; DCM.b(3,2,4)= 1; DCM.b(2,3,4)= 1;
DCM.c = [0 0 0 0 0; 1 1 0 1 0; 0 0 0 0 0];
DCM.d = zeros(3,3,0);

save(fullfile(save_path,'DCM_model_1_1x1.mat'),'DCM');






clear DCM
end