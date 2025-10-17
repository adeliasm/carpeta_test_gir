%% Script para encontrar muestra dentro de mi training y test del TFM sin dif signif de edad.

clc
clear
close all

% 
load ( '../results/Subset_Creation_Results/DEFSubsetTraining_17032023.mat', 'idx_trainp') 
load ( '../results/Subset_Creation_Results/DEFSubsetTest_17032023.mat', 'idx_testp')
load('D:\Postanoxic connectivity measures\ADELIA\Project_DCL\results\Subset_Creation_Results\TFM_finaldata_DMN_DCLvsControls_TFM.mat')

% load db original 
subj_datas = readtable('../data/DCLmegtusalen_v2.xlsx');
        % Diagnósticos: 
                    % 1 = Control (de MCI)
                    % 2 = control con QSM (QSM=queja subjetiva memoria)
                    % 3 = DCLa (a=amnésico)  % OJO: los que más riesgo tienen de evolucionar a AD
                    % 4 = DCLm (m=multi)
                    % 5 = DCLu (u=único)
                    % 6 = AD
                    % 7 = control con antecedentes de AD % proyecto FAM (tras calibrado)
                    % 8 = control sin QSM
                    % 9 = control sin antecedentes de AD % proyecto FAM (tras calibrado)
                    % 10 = No definido
subj_datas = subj_datas([find(subj_datas.diag == 1 | subj_datas.diag == 8 | subj_datas.diag == 3 | subj_datas.diag == 4 | subj_datas.diag == 5)],:);
subj_datas.Properties.VariableNames{1} = 'IdMEG';
subj_datas.IdMEG = strrep(subj_datas.IdMEG,'U1','UMEC-');
subj_datas.IdMEG = strrep(subj_datas.IdMEG,'N1','NEMOS-');
subj_datas = sortrows(subj_datas,'IdMEG','ascend'); %reordenar subj_datas para que tenga el mismo orden que plvdatas 
subj_datas(subj_datas.spectra_quality == 4, :) = []; % eliminar los que tienen registro de MEG malo (26 sujetos)
for i = 1:height(subj_datas)
     if  subj_datas.diag(i) == 1 || subj_datas.diag(i) == 8
         subj_datas.diag(i) = 0;
     else 
         subj_datas.diag(i) = 1;
     end 
end
db = subj_datas;
subj_datas = subj_datas(:,[1,4,6]);  % me quedo solo con ID, diag y age

        % comprobar train y test 
        Y = finaldata_DMN(1).table.outcome;
        nTrain = sum(idx_trainp) % 210 sujetos en training
        nTest = sum(idx_testp)  % 52 sujetos en test
        
        % saber proporción DCL en subsets training y test
             propDCL_train = sum(Y(idx_trainp))/length(Y(idx_trainp))*100  %55.24  (17/03/2023)
             proDCL_test  = sum(Y(idx_testp))/length(Y(idx_testp))*100      %55.77
                   


%crear subsamples train y test
subtrain = subj_datas(idx_trainp,:); 
subtest = subj_datas(idx_testp,:); 

train_controls = subtrain([find(subtrain.diag == 0)],:);
train_dcl = subtrain([find(subtrain.diag == 1)],:);

test_controls = subtest([find(subtest.diag == 0)],:);
test_dcl = subtest([find(subtest.diag == 1)],:);

%% seleccionar submuestra en train y test sin dif signif 

% ver dif edad: 
mean_ageControl = mean(train_controls.age) % = 70.49 // 70.51
std_ageControl = std(train_controls.age) % = 4.43 // 4.44
summary_train_controls = summary(train_controls(:,3))

mean_agetrain_dcl= mean(train_dcl.age) % = 74.06 // 74.12
std_agetrain_dcl = std(train_dcl.age) % = 5.12  // 4.97
summary_train_dcl = summary(train_dcl(:,3))  % age

figure;
subplot(1,2,1)
histogram(train_controls.age)
title('Edades Controles')

subplot(1,2,2)
histogram(train_dcl.age)
title('Edades DCL')

% test normalidad: 
addpath '../functions/swtest'
[H, pval, w] = swtest(train_controls.age)  % H = 0 --> distribución normal
[H, pval, w] = swtest(train_dcl.age) % H = 0 
[hT,pT,ciT,statsT] = ttest2(train_controls.age, train_dcl.age)

% % % CREAR SUBMUESTRA ALEATORIA PARA QUE NO HAYA DIFS SIGNIFS EN EDAD: 
% % % Parámetros iniciales
% % minAgeLimit = 64;
% % maxAgeLimit = 82;
% % minDur = 1;           % Mínima duración del rango de edad (por ejemplo 5 años)
% % maxIter = 5000;
% % 
% % % Inicialización
% % DCL_sample = [];
% % rng('default');       % Para reproducibilidad (puedes quitar si no lo necesitas)
% % 
% % % Bucle por rangos de edad
% % found = false;
% % 
% % for minAge = minAgeLimit:maxAgeLimit-1
% %     for maxAge = maxAgeLimit:-1:(minAge + minDur - 1)
% % 
% %         % Filtrar ambos grupos al rango actual
% %         Con2 = train_controls(train_controls.age >= minAge & train_controls.age <= maxAge, :);
% %         D2   = train_dcl   (train_dcl.age    >= minAge & train_dcl.age    <= maxAge, :);
% % 
% %         nC = height(Con2);
% %         nD = height(D2);
% % 
% %         % Saltar si no hay suficientes datos
% %         if nC < 5 || nD < nC
% %             continue
% %         end
% % 
% %         % Intentar submuestreo aleatorio
% %         for it = 1:maxIter
% %             sel = randperm(nD, nC);
% %             sampleD = D2(sel,:);
% %             % [h, p] = ttest2(Con2.age, sampleD.age);
% %             [p, h] = ranksum(Con2.age, sampleD.age);
% %             if h == 0
% %                 fprintf('✅ EXITO en rango [%d, %d] en iter %d, p=%.3f\n', ...
% %                         minAge, maxAge, it, p);
% %                 DCL_sample = sampleD;
% %                 found = true;
% %                 break;
% %             end
% %         end
% % 
% %         if found
% %             break;  % salir del segundo bucle
% %         end
% %     end
% %     if found
% %         break;  % salir del primer bucle
% %     end
% % end
% % 
% % % Resultado final
% % if ~found
% %     warning('❌ No se logró p>0.05 tras probar todos los rangos posibles.');
% % else
% %     % Puedes usar DCL_sample y Con2 para continuar tu análisis.
% %     fprintf('🎯 Submuestreo exitoso. N = %d\n', height(DCL_sample));
% % end
% % 
% % %%%---- Comprobación: 
% % ConAge = find(train_controls.age >= 65 & train_controls.age <= 79);
% % subsample_train_controls = train_controls(ConAge,:);
% % subsample_train_DCLs = train_dcl(ConAge,:);
% % [p, h] = ranksum(subsample_train_controls.age, subsample_train_DCLs.age)
% % 
% % % % ver si sirve tb para test: 
% % ConAge = find(test_controls.age >= 68 & test_controls.age <= 77);
% % subsample_test_controls = test_controls(ConAge,:);
% % subsample_test_DCLs = test_dcl(ConAge,:);
% % [p, h] = ranksum(subsample_test_controls.age, subsample_test_DCLs.age)  % NO SIRVE :(
% % 
% % % Buscar para test: 
% % minAgeLimit = 65;
% % maxAgeLimit = 79;
% % minDur = 1;           % Mínima duración del rango de edad (por ejemplo 5 años)
% % maxIter = 5000;
% % 
% % % Inicialización
% % DCL_sample = [];
% % rng('default');       % Para reproducibilidad (puedes quitar si no lo necesitas)
% % 
% % % Bucle por rangos de edad
% % found = false;
% % 
% % for minAge = minAgeLimit:maxAgeLimit-1
% %     for maxAge = maxAgeLimit:-1:(minAge + minDur - 1)
% % 
% %         % Filtrar ambos grupos al rango actual
% %         Con2 = test_controls(test_controls.age >= minAge & test_controls.age <= maxAge, :);
% %         D2   = test_dcl   (test_dcl.age    >= minAge & test_dcl.age    <= maxAge, :);
% % 
% %         nC = height(Con2);
% %         nD = height(D2);
% % 
% %         % Saltar si no hay suficientes datos
% %         if nC < 5 || nD < nC
% %             continue
% %         end
% % 
% %         % Intentar submuestreo aleatorio
% %         for it = 1:maxIter
% %             sel = randperm(nD, nC);
% %             sampleD = D2(sel,:);
% %             % [h, p] = ttest2(Con2.age, sampleD.age);
% %             [p, h] = ranksum(Con2.age, sampleD.age);
% %             if h == 0
% %                 fprintf('✅ EXITO en rango [%d, %d] en iter %d, p=%.3f\n', ...
% %                         minAge, maxAge, it, p);
% %                 DCL_sample = sampleD;
% %                 found = true;
% %                 break;
% %             end
% %         end
% % 
% %         if found
% %             break;  % salir del segundo bucle
% %         end
% %     end
% %     if found
% %         break;  % salir del primer bucle
% %     end
% % end
% % 
% % % Resultado final
% % if ~found
% %     warning('❌ No se logró p>0.05 tras probar todos los rangos posibles.');
% % else
% %     % Puedes usar DCL_sample y Con2 para continuar tu análisis.
% %     fprintf('🎯 Submuestreo exitoso. N = %d\n', height(DCL_sample));
% % end
% % 
% % 
% % %%%---- Comprobación: 
% % ConAge = find(test_controls.age >= 64 & test_controls.age <= 82);
% % subsample_test_controls = test_controls(ConAge,:);
% % subsample_test_DCLs = DCL_sample; 
% % [p, h] = ranksum(subsample_test_controls.age, subsample_test_DCLs.age)
% % 
% % % ver si hay dif entre train y test: 
% % [p, h] = ranksum(subsample_test_controls.age, subsample_train_controls.age)  % NO SIRVE :(
% % [p, h] = ranksum(subsample_test_DCLs.age, subsample_train_DCLs.age)  



%% new try: 
addpath '../functions/'; 

minAgeLimit = 64;
maxAgeLimit = 82;
minDur = 5;
maxIter = 5000;

% --- Train ---
[train_controls_sub, train_dcl_sub] = match_age_groups(train_controls, train_dcl, minAgeLimit, maxAgeLimit, minDur, maxIter);
    % solución: 68-77
% --- Test ---
[test_controls_sub, test_dcl_sub] = match_age_groups(test_controls, test_dcl, minAgeLimit, maxAgeLimit, minDur, maxIter);
    % solución: 64-82

fprintf('\nVerificando diferencias de edad...\n');

% 1. Entre controles y DCL en TRAIN
[p1, h1] = ranksum(train_controls_sub.age, train_dcl_sub.age);
fprintf('TRAIN: Controles vs DCL, p=%.4f\n', p1);

% 2. Entre controles y DCL en TEST
[p2, h2] = ranksum(test_controls_sub.age, test_dcl_sub.age);
fprintf('TEST: Controles vs DCL, p=%.4f\n', p2);

% 3. Entre TRAIN y TEST en controles
[p3, h3] = ranksum(train_controls_sub.age, test_controls_sub.age);
fprintf('TRAIN vs TEST: Controles, p=%.4f\n', p3);

% 4. Entre TRAIN y TEST en DCL
[p4, h4] = ranksum(train_dcl_sub.age, test_dcl_sub.age);
fprintf('TRAIN vs TEST: DCL, p=%.4f\n', p4);

summary(train_controls_sub(:,3))
summary(train_dcl_sub(:,3))
summary(test_controls_sub(:,3))
summary(test_dcl_sub(:,3))

% aleatorizar las filas
train_age_check = [train_controls_sub; train_dcl_sub]; 
rand = randperm(height(train_age_check)); 
train_age_check = train_age_check(rand, :);  
test_age_check = [test_controls_sub; test_dcl_sub]; 
rand = randperm(height(test_age_check)); 
test_age_check = test_age_check(rand, :);  

% Guardar nuevo FD con age_checked
% new age_checked: 
load('D:\Postanoxic connectivity measures\ADELIA\Project_DCL\results\Subset_Creation_Results\TFM_finaldata_DMN_DCLvsControls_TFM.mat')
load('D:\Postanoxic connectivity measures\ADELIA\Project_DCL\derivatives\TFM_test_age_check.mat')
load('D:\Postanoxic connectivity measures\ADELIA\Project_DCL\derivatives\TFM_train_age_check.mat')

% %modif para usar las mismas 45 feats que en mi TFM. 
% load('D:\Postanoxic connectivity measures\ADELIA\Project_DCL\results\Features_Lasso1_TFM_beta.mat')
% % load('D:\Postanoxic connectivity measures\ADELIA\Project_DCL\results\dataClassif_TFM_beta.mat')
% 
% feats(1) = []; % elimino la columna 'Group'
% [~, ind_feat_final] = ismember(feats, finaldata_DMN(8).table.Properties.VariableNames);  % saca los idx de las cols de finaldata manteniendo el orden de 'feats'
% 
% data_ = [table(finaldata_DMN(8).table.IdMEG, finaldata_DMN(8).table.outcome, 'VariableNames', {'IdMEG', 'Group'}), ...
%     finaldata_DMN(8).table(:, ind_feat_final)];

% seleccionar las filas de los sujs que quiero
checked_sujs = [train_age_check.IdMEG; test_age_check.IdMEG];
finaldata_DMN_0 = finaldata_DMN;

for bindex = 1: numel(finaldata_DMN)
   
    [~, idx_sujs] = ismember(checked_sujs, finaldata_DMN(bindex).table.IdMEG);
    finaldata_DMN(bindex).table = finaldata_DMN(bindex).table(idx_sujs,:);

    [idx_trainp, ~] = ismember(finaldata_DMN(bindex).table.IdMEG, train_age_check.IdMEG);  % 136 sujs
    [idx_testp, ~] = ismember(finaldata_DMN(bindex).table.IdMEG, test_age_check.IdMEG);   % 44 sujs

    nTrain = sum(idx_trainp) % 136 sujetos en training
    nTest = sum(idx_testp)  % 44 sujetos en test

    % saber proporción DCL en subsets training y test
    Y = finaldata_DMN(bindex).table.outcome;
    prop1_train = sum(Y(idx_trainp))/length(Y(idx_trainp))*100  %50%
    prop1_test  = sum(Y(idx_testp))/length(Y(idx_testp))*100      %50%
end

save('../results/Subset_Creation_Results/FD_TFM_age_checked_180sujs_DCLvsControls_30072025.mat', 'finaldata_DMN')
save('../results/Subset_Creation_Results/SubsetTraining_TFM_age_checked_180sujs_DCLvsControls_30072025.mat', 'idx_trainp')
save('../results/Subset_Creation_Results/SubsetTest_TFM_age_checked_180sujs_DCLvsControls_30072025.mat', 'idx_testp')