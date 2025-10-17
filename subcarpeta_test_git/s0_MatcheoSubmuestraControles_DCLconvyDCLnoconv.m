%% EXPLORADOR DE DATOS DCL
clc
clear 
close all

db_megtusalen = readtable('../data/DCLmegtusalen_v2.xlsx');
db_DCLconv = readtable('../../ICON/derivatives/22052025_dbDCLconv.xlsx'); 
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


db_megtusalen = db_megtusalen([find(db_megtusalen.diag == 1 | db_megtusalen.diag == 8 )],:);
db_megtusalen = sortrows(db_megtusalen,'Code','ascend'); %reordenar subj_datas para que tenga el mismo orden que plvdatas 
idx_umec = startsWith(db_megtusalen.Code, 'U', 'IgnoreCase', true);
db_megtusalen(idx_umec,:)= []; % eliminar los que son de umec, me quedo solo con los de nemos. 
db_megtusalen(db_megtusalen.spectra_quality == 4, :) = []; % eliminar los que tienen registro de MEG malo (26 sujetos)

for i = 1:height(db_megtusalen)
     if  db_megtusalen.diag(i) == 1 || db_megtusalen.diag(i) == 8
         db_megtusalen.diag(i) = 0;
     else 
         db_megtusalen.diag(i) = 1;
     end 
end

controls = db_megtusalen; 
pDCL = db_DCLconv(db_DCLconv.conv == 1, :);
sDCL = db_DCLconv(db_DCLconv.conv == 0, :);

%% VOY POR AQUÍ, ME QUEDA SEGUIR MODIFICANDO EL SCRIPT, LO ÚLTIMO QUE TENGO ES LA SEPRACIÓN DE CONTROLES, pDCL y sDCL. --- 29092025

% histogram(subj_datas.diag)
% title('Controles vS pDCL');
% xlabel('');
% ylabel('');
% 
% % % Comprobar diagnósticos de controles y pDCL de ambos proyectos según MMSE
% % 
% % a = find(controls.mmse < 27); % Hay uno con mmse = 26, pero hay que tener en cuenta edad (=71) y años de educación (=9), por lo que según la tabla de mmse corregido para edad y años de educación seguiría siendo sano. 
% % b = find(pDCL.mmse > 27); % Hay 56 sujetos clasificados como pDCL que tienen mmse > 27, PERO EN OTRAS PRUEBAS DE MEMORIA TIENEN PUNTUCACIONES BAJISIMAS !!
% % c = pDCL(b,:); 
% % 
% % a = find(isnan(controls.mmse)); % No hay NaN de mmse en controles
% % b = find(isnan(pDCL.mmse)); % Hay 6 NaN de mmse en pDCL
% % 
% % a = find(isnan(controls.age)); % No hay NaN 
% % b = find(isnan(pDCL.age));
% % 
% % % % % Replace mmse NaN by median value from all subjects
% % % % pDCL.mmse(isnan(pDCL.mmse)) = median(subj_datas.mmse, 'omitnan');
% % % % 
% % 
% % find(isnan(controls.edu_years)) % Hay 1 NaN de edu_years en controls
% % find(isnan(pDCL.edu_years)) % Hay 9 NaN de edu_years en DCL
% % 
% % 
% % % Replace edu_years NaN by median value from all subjects
% % controls.edu_years(isnan(controls.edu_years)) = median(subj_datas.edu_years, 'omitnan');
% % pDCL.edu_years(isnan(pDCL.edu_years)) = median(subj_datas.edu_years, 'omitnan');
% 
% %% ELIMINAR LOS IS NAN
% a = find(isnan(subj_datas.mmse));
% subj_datas([a], :) = []; 
% a = find(isnan(subj_datas.edu_years)); 
% subj_datas([a], :) = []; % quedan 132 pDCL y 116 controles
% 
% controls =  subj_datas([find(subj_datas.diag == 0)],:); %116
% pDCL = subj_datas([find(subj_datas.diag ==1)],:);  %132
% 
% %% // PROBAR A ELIMINAR LOS NAN
% 
% 
% %% Tests estadísticos 
% 
% % Saber si siguen distribución normal
%     %  Shapiro-Wilk parametric hypothesis test of composite normality
%     % The Shapiro-Wilk and Shapiro-Francia null hypothesis is: 
%     %   "X is normal with unspecified mean and variance."
% addpath '../functions/swtest'
% % [H, pval, w] = swtest(controls.mmse)  % H = 1 --> distribución no normal
% % [H, pval, w] = swtest(pDCL.mmse) % H = 1 --> distribución no normal
% [H, pval, w] = swtest(subj_datas.mmse) % H = 1  --> distribución no normal
% 
% % Saber si hay diferencias significativas de MMSE entre los 2 grupos
%     % Mann-Whitney-Wilcoxon non parametric test for two unpaired groups.
% % addpath '../functions/mwwtest'
% % mwwtest(controls.mmse, pDCL.mmse)  % H=1 --> signif, pval = 0.00000 (< 0.001)
% % [p,h] = ranksum(controls.mmse, pDCL.mmse) 
% 
% 
% 
% aoctool (subj_datas.mmse, subj_datas.age, subj_datas.diag)   % F = 4.1, pval = 0.044 // F=4.17, p= 0.0423
% 
% % El análisis de covarianza permite aumentar la precisión de los experimentos y eliminar los 
% % efectos de variables que no tienen nada que ver con el tratamiento, pero que sin embargo, sí están influyendo en los resultados.
% 
% mean_mmseControl = mean(controls.mmse) % = 29 // 28.99
% std_mmseControl = std(controls.mmse) % = 1.1 // 1.1
% mean_mmseDCL = mean(pDCL.mmse) % = 26.66 // 26.63
% std_mmseDCL = std(pDCL.mmse) % = 2.52 // 2.6
% 
% histogram(controls.mmse)
% title('Puntuaciones MMSE');
% xlabel('puntuación');
% ylabel(''); 
% 
% hold on
% 
% histogram(pDCL.mmse)
% legend ('controls', 'pDCL')
% 
% histogram(subj_datas.mmse)
% 
% hold off
% 
   
% Comprobar diferencias EDAD entre ambos grupos

% ageControlNEMOS = mean(controls.age(1:56));    % 70.20
% ageControlUMEC = mean(controls.age(56:end));   % 70.85
mean_ageControl = mean(controls.age) % = 70.49 // 70.51
std_ageControl = std(controls.age) % = 4.43 // 4.44
summary_controls = summary(controls(:,6))  % age

% ageDCLNEMOS = mean(pDCL.age(1:56));     % 74.75
% ageDCLUMEC = mean(pDCL.age(56:end));    % 73.68
mean_agepDCL = mean(pDCL.age) % = 74.06 // 74.12
std_agepDCL = std(pDCL.age) % = 5.12  // 4.97
summary_pDCL = summary(pDCL(:,7))  % age

mean_agesDCL= mean(sDCL.age) % = 74.06 // 74.12
std_agesDCL = std(sDCL.age) % = 5.12  // 4.97
summary_sDCL = summary(sDCL(:,7))  % age

figure;
subplot(1,2,1)
histogram(controls_sample.age)
title('Edades Controles')

subplot(1,2,2)
histogram(DCL_sample.age)
title('Edades pDCL')

subplot(1,3,3)
histogram(sDCL.age)
title('Edades pDCL')

% Saber si siguen distribución normal
    %  Shapiro-Wilk parametric hypothesis test of composite normality
    % The Shapiro-Wilk and Shapiro-Francia null hypothesis is: 
    %   "X is normal with unspecified mean and variance."
addpath '../functions/swtest'
[H, pval, w] = swtest(controls.age)  % H = 0 --> distribución normal
[H, pval, w] = swtest(pDCL.age) % H = 0 
% [H, pval, w] = swtest(subj_datas.age)  % H = 0 
[hT,pT,ciT,statsT] = ttest2(controls.age, sDCL.age)  % H = 1, pval =  8.4227e-09 (<0.001)
% mwwtest(controls.age, pDCL.age) % H = 1


% CREAR SUBMUESTRA ALEATORIA PARA QUE NO HAYA DIFS SIGNIFS EN EDAD: 
% Parámetros iniciales
minAgeLimit = 60;
maxAgeLimit = 87;
minDur = 1;           % Mínima duración del rango de edad (por ejemplo 5 años)
maxIter = 5000;

% Inicialización
pDCL_sample = [];
rng('default');       % Para reproducibilidad (puedes quitar si no lo necesitas)

% Bucle por rangos de edad
found = false;

for minAge = minAgeLimit:maxAgeLimit-1
    for maxAge = maxAgeLimit:-1:(minAge + minDur - 1)
        
        % Filtrar ambos grupos al rango actual
        Con2 = controls(controls.age >= minAge & controls.age <= maxAge, :);
        D2   = pDCL   (pDCL.age    >= minAge & pDCL.age    <= maxAge, :);
        
        nC = height(Con2);
        nD = height(D2);
        
        % Saltar si no hay suficientes datos
        if nC < 5 || nD < nC
            continue
        end
        
        % Intentar submuestreo aleatorio
        for it = 1:maxIter
            sel = randperm(nD, nC);
            sampleD = D2(sel,:);
            % [h, p] = ttest2(Con2.age, sampleD.age);
            [p, h] = ranksum(Con2.age, sampleD.age);
            if h == 0
                fprintf('✅ EXITO en rango [%d, %d] en iter %d, p=%.3f\n', ...
                        minAge, maxAge, it, p);
                pDCL_sample = sampleD;
                found = true;
                break;
            end
        end
        
        if found
            break;  % salir del segundo bucle
        end
    end
    if found
        break;  % salir del primer bucle
    end
end

% Resultado final
if ~found
    warning('❌ No se logró p>0.05 tras probar todos los rangos posibles.');
else
    % Puedes usar pDCL_sample y Con2 para continuar tu análisis.
    fprintf('🎯 Submuestreo exitoso. N = %d\n', height(pDCL_sample));
end

%% Modificación chatgpt para N > 38 sujs: 
% --- Parámetros iniciales ---
minAgeLimit = 60;
maxAgeLimit = 87;
minDur = 1;           % duración mínima del rango
maxIter = 5000;       % nº máximo de intentos por rango
nMinControls = 32;    % tamaño mínimo de controles

% --- Inicialización ---
controls_sample = [];
rng('default');       % Para reproducibilidad

% --- Bucle por rangos de edad ---
found = false;

for minAge = minAgeLimit:maxAgeLimit-1
    for maxAge = maxAgeLimit:-1:(minAge + minDur - 1)
        
        % Filtrar ambos grupos al rango actual
        Con2 = controls(controls.age >= minAge & controls.age <= maxAge, :);
        D2   = sDCL(sDCL.age >= minAge & sDCL.age <= maxAge, :);
        
        nC = height(Con2); % controles en rango
        nD = height(D2);   % DCL en rango
        
        % Necesitamos al menos 38 controles y los 41 DCL (o los que queden)
        if nC < nMinControls || nD < 30 % por si hay pocos DCL
            continue
        end
        
        % Intentar submuestreo aleatorio
        for it = 1:maxIter
            
            % Selecciona al azar controles (al menos 38)
            nSel = nMinControls;  % puedes cambiarlo a nC si quieres usar todos
            selC = randperm(nC, nSel);
            sampleC = Con2(selC,:);
            
            % Aquí comparas la edad entre sDCL (o D2) y los controles elegidos
            % [p, h] = ranksum(D2.age, sampleC.age); % p>0.05 = sin dif
            [h,p] = ttest2(D2.age, sampleC.age);  % alternativa paramétrica
            
            if h == 0 && p > 0.05
                fprintf('✅ EXITO rango [%d, %d] iter %d, p=%.3f\n', ...
                        minAge, maxAge, it, p);
                controls_sample = sampleC;
                DCL_sample = D2;   % por si quieres guardar también DCL filtrados
                found = true;
                break;
            end
        end
        
        if found
            break;  % salir del segundo bucle
        end
    end
    if found
        break;  % salir del primer bucle
    end
end

% --- Resultado final ---
if ~found
    warning('❌ No se logró p>0.05 tras probar todos los rangos posibles.');
else
    fprintf('🎯 Submuestreo exitoso. Controles N=%d, DCL N=%d\n', ...
            height(controls_sample), height(DCL_sample));
end

%% Ahora matchear también con controles y pDCLs que ya seleccioné: 
db = readtable('D:\PROGRAMMING\MATLAB\Postanoxic connectivity measures\ADELIA\Project_DCL\derivatives\AdeliaDB_pDCLvsControls_agebalanced_29092025.xlsx');
controls = db(db.conv == 0,:); 
pDCL = db(db.conv ==1, :); 

% --- Parámetros iniciales ---
minAgeLimit = 60;
maxAgeLimit = 87;
minDur = 1;           % duración mínima del rango
maxIter = 5000;       % nº máximo de intentos por rango
nSel = 32;            % tamaño deseado de la submuestra sDCL

% --- Inicialización ---
rng('default');       % Para reproducibilidad
found = false;

% --- Bucle por rangos de edad ---
for minAge = minAgeLimit:maxAgeLimit-1
    for maxAge = maxAgeLimit:-1:(minAge + minDur - 1)
        
        % Filtrar grupos al rango actual
        Con2 = controls(controls.age >= minAge & controls.age <= maxAge, :);
        pD2  = pDCL(pDCL.age >= minAge & pDCL.age <= maxAge, :);
        sD2  = sDCL(sDCL.age >= minAge & sDCL.age <= maxAge, :);
        
        nC = height(Con2);
        nP = height(pD2);
        nS = height(sD2);
        
        % Necesitamos al menos 32 en cada grupo
        if nC < nSel || nP < nSel || nS < nSel
            continue
        end
        
        % Intentar submuestreo aleatorio de sDCL
        for it = 1:maxIter
            
            selS = randperm(nS, nSel);      % submuestra de sDCL
            sampleS = sD2(selS,:);
            
            % Comprobamos dos tests de edad:
            [~,p1] = ttest2(sampleS.age, Con2.age); % sDCL vs controles
            [~,p2] = ttest2(sampleS.age, pD2.age);  % sDCL vs pDCL
            
            % O si prefieres no paramétrico:
            % p1 = ranksum(sampleS.age, Con2.age);
            % p2 = ranksum(sampleS.age, pD2.age);
            
            if p1 > 0.05 && p2 > 0.05
                fprintf('✅ EXITO rango [%d, %d] iter %d, p1=%.3f, p2=%.3f\n', ...
                        minAge, maxAge, it, p1, p2);
                sDCL_sample = sampleS;
                controls_sample = Con2;
                pDCL_sample = pD2;
                found = true;
                break;
            end
        end
        
        if found, break; end
    end
    if found, break; end
end

% --- Resultado final ---
if ~found
    warning('❌ No se logró p>0.05 tras probar todos los rangos posibles.');
else
    fprintf('🎯 Submuestreo exitoso. Controles=%d, pDCL=%d, sDCL=%d\n', ...
            height(controls_sample), height(pDCL_sample), height(sDCL_sample));
end


%% ---- Comprobación: 
ConAge = find(controls.age >= 69 & controls.age <= 80);
subsample_Controls = controls(ConAge,:);
subsample_DCLs = pDCL_sample; 
[hT,pT] = ttest2(controls_sample.age, sDCL_sample.age)  % h = 0, NO SIGNIF !!
[hT,pT] = ttest2(pDCL_sample.age, sDCL_sample.age)  % h = 0, NO SIGNIF !!

mean_ageControl = mean(subsample_Controls.age) % = 70.49 // 70.51
std_ageControl = std(subsample_Controls.age) % = 4.43 // 4.44
summary_controls = summary(subsample_Controls(:,6))  % age

mean_ageDCLs= mean(subsample_DCLs.age) % = 74.06 // 74.12
std_ageDCLs = std(subsample_DCLs.age) % = 5.12  // 4.97
summary_dcls = summary(subsample_DCLs(:,6))  % age

figure;
subplot(1,2,1)
histogram(subsample_Controls.age)
title('Edades subsample Controles')

subplot(1,2,2)
histogram(subsample_DCLs.age)
title('Edades subsample DCL')

addpath '../functions/swtest'
[H, pval, w] = swtest(subsample_Controls.age)  % H = 1 --> distribución NO normal
[H, pval, w] = swtest(subsample_DCLs.age) % H = 1 
[p, h] = ranksum(subsample_Controls.age, subsample_DCLs.age)

pDCL_sample = DCL_sample;
cols = {'Code', 'conv', 'spectra_quality', 'conv_time', 'dx', 'age', 'sex', 'years_edu', 'Nivel_Escolaridad_v2', 'MMSE'}; 
idx_cols = ismember(DCL_sample.Properties.VariableNames, cols);
pDCL_sample = pDCL_sample(:, idx_cols);

cols = {'Code', 'conv', 'spectra_quality', 'conv_time', 'dx', 'age', 'sex', 'years_edu', 'Nivel_Escolaridad_v2', 'MMSE'}; 
idx_cols = ismember(sDCL_sample.Properties.VariableNames, cols);
sDCL_sample = sDCL_sample(:, idx_cols);

cols = {'Code', 'conv', 'spectra_quality', 'diag', 'age', 'sex', 'edu_years', 'edu_clas_level', 'mmse'}; 
idx_cols = ismember(controls_sample.Properties.VariableNames, cols); 
controls_sample = controls_sample(:,idx_cols); 

sDCL_sample.Properties.VariableNames =  {'Code', 'conv', 'spectra_quality', 'conv_time', 'diag', 'age', 'sex', 'edu_years', 'edu_clas_level', 'mmse'}; 
controls_sample.conv_time = zeros(height(controls_sample),1);

finaldb = cat(1, controls_sample, pDCL_sample);
finaldb = cat(1,finaldb, sDCL_sample); 
writetable(finaldb, '../derivatives/AdeliaDB_sDCLvs_pDCLvsControls_agebalanced_29092025.xlsx')



%%%%%% 
% [hC,pC,ciC,statsC] = ttest2(controls.age(1:56),controls.age(56:end))
% [hD,pD,ciD,statsD] = ttest2(pDCL.age(1:56),pDCL.age(56:end))               % = 0 por lo que no hay dif estadísticas de edad de controles ni DCL ENTRE AMBOS PROYECTOS
% [hT,pT,ciT,statsT] = ttest2(controls.age,pDCL.age)              % = 1, por lo que SÍ HAY dif estadísticamente significativas entre edad de CONTROLES VS pDCL en general. 
%     % The result h is 1 if the test rejects the null hypothesis at the 5% significance level, and 0 otherwise.
%     % Test the null hypothesis that the two data samples are from populations with equal means.
% histogram(subj_datas.age)% 

% histogram(controls.age)
% title('Edad');
% xlabel('años');
% ylabel(''); 
% 
% hold on
% 
% histogram(pDCL.age)
% legend ('controls', 'pDCL')
% 
% hold off
% 
% 
% x = [0 1];
% data = [mean_ageDCLs mean_ageControl];
% errhigh = std_ageControl;
% errlow  = [4.4 2.4 2.3 0.5 1.6 1.5 4.5 1.5 0.4 1.2 1.3 0.8 1.9];
% 
% bar(x,data)                
% 
% hold on
% 
% bar(mean_ageDCLs, mean_ageControl)
% hold on
% er = errorbar(x,std_ageControl,errlow,errhigh);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% 
% hold off


% Explorar distribución por SEXO
% 1 = hombre, 2 = mujer 
% female = subj_datas(find(subj_datas.sex == 2), :); % 168 mujeres --> 92 pDCL, 76 controles // 157: 82 dcl, 75 controles
% male = subj_datas(find(subj_datas.sex == 1), :); % 94 hombres --> 53 pDCL, 41 controles // 91: 50 dcl, 41 controles
% sum(female.diag)
% sum(male.diag)

% Fisher's exact test for gender (paper sandra pusil)
       % info: fishertest(x) returns a test decision for Fisheris exact test of the null hypothesis that 
        % there are no nonrandom associations between the two categorical variables in x, against 
        % the alternative that there is a nonrandom association. The result h is 1 if the test rejects the 
        % null hypothesis at the 5% significance level, or 0 otherwise.
        % FISHERTEST requires a 2-by-2 matrix or table X.
% 
% sex = table([75;41], [82;50], 'VariableNames',{'Control','DCL'},'RowNames',{'Female','Male'})    
% sex = table([76;41], [92;53], 'VariableNames',{'Control','DCL'},'RowNames',{'Female','Male'}) % p = 0.897

sex = table([table2array(sum(finaldb(1:76,7) == 2)); table2array(sum(finaldb(1:76,7) == 1))],...
    [table2array(sum(finaldb(77:end,7) == 2)); table2array(sum(finaldb(77:end,7) == 1))],...
    'VariableNames',{'DCL','Control'},'RowNames',{'Female','Male'}) % p = 0.897
[h,p,stats] = fishertest(sex) 
            % H= 0 --> No hay diferencias significativas en la distribución por sexo entre tus grupos DCL y Control.

% 
% % Explorar distribución por AÑOS DE EDUCACIÓN
% histogram(subj_datas.edu_years)
% histogram(controls.edu_years)
% title('Años de educación');
% xlabel('años');
% ylabel(''); 
% 
% hold on
% 
% histogram(pDCL.edu_years)
% legend ('controls', 'pDCL')
% 
% hold off
% % 
% % [H, pval, w] = swtest(controls.edu_years)  % H = 1 --> distribución no normal
% % [H, pval, w] = swtest(pDCL.edu_years) % H = 1 
% [H, pval, w] = swtest(subj_datas.edu_years) % H = 1  --> distribución no normal
% 
% % Saber si hay diferencias significativas entre los 2 grupos
%     % Mann-Whitney-Wilcoxon non parametric test for two unpaired groups.
% 
% % mwwtest(controls.edu_years, pDCL.edu_years)  % H=1 --> signif, pval = (<0.001)
% 
% aoctool (subj_datas.edu_years, subj_datas.age, subj_datas.diag)   % F = 0.02, pval = 0.885 //F=0.01, p=0.937
% 
% mean_eduControl = mean(controls.edu_years) % = 13.65
% std_eduControl = std(controls.edu_years) % = 5.24
% mean_eduDCL = mean(pDCL.edu_years) % = 8.98
% std_eduDCL = std(pDCL.edu_years) % = 4.37
% 
% 
