classdef IndexConst
  properties (Constant)
    % Indices for results arrays, Xinfo, Yinfo and ModelInfo
    NUM_VARS = 3;
    NOISE_IND = 1;
    CLEAN_IND = 2;
    RECON_IND = 3;
    B_IND = 1;
    PHI_IND = 2;
    W_IND = 3;
    
    % Indices for runtime variables
    RuntimeNames = {'Training','Recovery'};
    TRAINING_TIME_IND = 1;
    RECOVERY_TIME_IND = 2;

    % Indices for objective values    
    ObjectiveValueNames = {'Objective', 'L1 Loss', 'L2 Loss', 'Reg Loss'}; 
    OBJ_IND = 1;
    L1_IND = 2;
    L2_IND = 3;
    REG_IND = 4;
  end
end