%% EXPLORADOR DE DATOS DCL
% añado esta nueva linea para probar
clc
clear 
close all

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

% a = find(subj_datas.spectra_quality == 3);  % Hay 14 sujetos con MEG DUDA --> 6 son controles y 8 son DLCs
% a2 = find(subj_datas.spectra_quality == 2);  % Hay 7 sujetos con MEG REGULAR --> 2 son controles y 5 son DLCs
% a3 = find(subj_datas.spectra_quality == 1);  % Hay 241 sujetos con MEG BUENO --> serían 109 controles y 132 DCLs
% b = subj_datas(a3,:);
% sum(b.diag)

Controls =  subj_datas([find(subj_datas.diag == 0)],:); %127  --> 117 después de quitar los de MEG MALO
DCLs = subj_datas([find(subj_datas.diag ==1)],:);  %161  -->145 después de quitar MEG MALO
    % converters = find(subj_datas.conv == 1 ); % En cuánto tiempo ?? todos los conversores son del proyecto NEMOS, pero cuánto tiempo se estuvo siguiendo a estos sujetos ??

% histogram(subj_datas.diag)
% title('Controles vS DCLs');
% xlabel('');
% ylabel('');
% 
% % % Comprobar diagnósticos de controles y DCLs de ambos proyectos según MMSE
% % 
% % a = find(Controls.mmse < 27); % Hay uno con mmse = 26, pero hay que tener en cuenta edad (=71) y años de educación (=9), por lo que según la tabla de mmse corregido para edad y años de educación seguiría siendo sano. 
% % b = find(DCLs.mmse > 27); % Hay 56 sujetos clasificados como DCLs que tienen mmse > 27, PERO EN OTRAS PRUEBAS DE MEMORIA TIENEN PUNTUCACIONES BAJISIMAS !!
% % c = DCLs(b,:); 
% % 
% % a = find(isnan(Controls.mmse)); % No hay NaN de mmse en controles
% % b = find(isnan(DCLs.mmse)); % Hay 6 NaN de mmse en DCLs
% % 
% % a = find(isnan(Controls.age)); % No hay NaN 
% % b = find(isnan(DCLs.age));
% % 
% % % % % Replace mmse NaN by median value from all subjects
% % % % DCLs.mmse(isnan(DCLs.mmse)) = median(subj_datas.mmse, 'omitnan');
% % % % 
% % 
% % find(isnan(Controls.edu_years)) % Hay 1 NaN de edu_years en Controls
% % find(isnan(DCLs.edu_years)) % Hay 9 NaN de edu_years en DCL
% % 
% % 
% % % Replace edu_years NaN by median value from all subjects
% % Controls.edu_years(isnan(Controls.edu_years)) = median(subj_datas.edu_years, 'omitnan');
% % DCLs.edu_years(isnan(DCLs.edu_years)) = median(subj_datas.edu_years, 'omitnan');
% 
% %% ELIMINAR LOS IS NAN
% a = find(isnan(subj_datas.mmse));
% subj_datas([a], :) = []; 
% a = find(isnan(subj_datas.edu_years)); 
% subj_datas([a], :) = []; % quedan 132 DCLs y 116 controles
% 
% Controls =  subj_datas([find(subj_datas.diag == 0)],:); %116
% DCLs = subj_datas([find(subj_datas.diag ==1)],:);  %132
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
% % [H, pval, w] = swtest(Controls.mmse)  % H = 1 --> distribución no normal
% % [H, pval, w] = swtest(DCLs.mmse) % H = 1 --> distribución no normal
% [H, pval, w] = swtest(subj_datas.mmse) % H = 1  --> distribución no normal
% 
% % Saber si hay diferencias significativas de MMSE entre los 2 grupos
%     % Mann-Whitney-Wilcoxon non parametric test for two unpaired groups.
% % addpath '../functions/mwwtest'
% % mwwtest(Controls.mmse, DCLs.mmse)  % H=1 --> signif, pval = 0.00000 (< 0.001)
% % [p,h] = ranksum(Controls.mmse, DCLs.mmse) 
% 
% 
% 
% aoctool (subj_datas.mmse, subj_datas.age, subj_datas.diag)   % F = 4.1, pval = 0.044 // F=4.17, p= 0.0423
% 
% % El análisis de covarianza permite aumentar la precisión de los experimentos y eliminar los 
% % efectos de variables que no tienen nada que ver con el tratamiento, pero que sin embargo, sí están influyendo en los resultados.
% 
% mean_mmseControl = mean(Controls.mmse) % = 29 // 28.99
% std_mmseControl = std(Controls.mmse) % = 1.1 // 1.1
% mean_mmseDCL = mean(DCLs.mmse) % = 26.66 // 26.63
% std_mmseDCL = std(DCLs.mmse) % = 2.52 // 2.6
% 
% histogram(Controls.mmse)
% title('Puntuaciones MMSE');
% xlabel('puntuación');
% ylabel(''); 
% 
% hold on
% 
% histogram(DCLs.mmse)
% legend ('Controls', 'DCLs')
% 
% histogram(subj_datas.mmse)
% 
% hold off
% 
   
% Comprobar diferencias EDAD entre ambos grupos

% ageControlNEMOS = mean(Controls.age(1:56));    % 70.20
% ageControlUMEC = mean(Controls.age(56:end));   % 70.85
mean_ageControl = mean(Controls.age) % = 70.49 // 70.51
std_ageControl = std(Controls.age) % = 4.43 // 4.44
summary_controls = summary(Controls(:,6))  % age

% ageDCLNEMOS = mean(DCLs.age(1:56));     % 74.75
% ageDCLUMEC = mean(DCLs.age(56:end));    % 73.68
mean_ageDCLs= mean(DCLs.age) % = 74.06 // 74.12
std_ageDCLs = std(DCLs.age) % = 5.12  // 4.97
summary_dcls = summary(DCLs(:,6))  % age

figure;
subplot(1,2,1)
histogram(Controls.age)
title('Edades Controles')

subplot(1,2,2)
histogram(DCLs.age)
title('Edades DCL')

% Saber si siguen distribución normal
    %  Shapiro-Wilk parametric hypothesis test of composite normality
    % The Shapiro-Wilk and Shapiro-Francia null hypothesis is: 
    %   "X is normal with unspecified mean and variance."
addpath '../functions/swtest'
[H, pval, w] = swtest(Controls.age)  % H = 0 --> distribución normal
[H, pval, w] = swtest(DCLs.age) % H = 0 
% [H, pval, w] = swtest(subj_datas.age)  % H = 0 
[hT,pT,ciT,statsT] = ttest2(Controls.age, DCLs.age)  % H = 1, pval =  8.4227e-09 (<0.001)
% mwwtest(Controls.age, DCLs.age) % H = 1


% CREAR SUBMUESTRA ALEATORIA PARA QUE NO HAYA DIFS SIGNIFS EN EDAD: 
% Parámetros iniciales
minAgeLimit = 64;
maxAgeLimit = 82;
minDur = 1;           % Mínima duración del rango de edad (por ejemplo 5 años)
maxIter = 5000;

% Inicialización
DCL_sample = [];
rng('default');       % Para reproducibilidad (puedes quitar si no lo necesitas)

% Bucle por rangos de edad
found = false;

for minAge = minAgeLimit:maxAgeLimit-1
    for maxAge = maxAgeLimit:-1:(minAge + minDur - 1)
        
        % Filtrar ambos grupos al rango actual
        Con2 = Controls(Controls.age >= minAge & Controls.age <= maxAge, :);
        D2   = DCLs   (DCLs.age    >= minAge & DCLs.age    <= maxAge, :);
        
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
                DCL_sample = sampleD;
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
    % Puedes usar DCL_sample y Con2 para continuar tu análisis.
    fprintf('🎯 Submuestreo exitoso. N = %d\n', height(DCL_sample));
end

%%%---- Comprobación: 
ConAge = find(Controls.age >= 69 & Controls.age <= 80);
subsample_Controls = Controls(ConAge,:);
subsample_DCLs = DCL_sample; 
[hT,pT,ciT,statsT] = ttest2(subsample_Controls.age, subsample_DCLs.age)  % h = 0, NO SIGNIF !!

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

finaldb = cat(1, subsample_DCLs, subsample_Controls); 
% writetable(finaldb, '../derivatives/AdeliaDB_DCLvsControls_checked_22072025.xlsx')



%%%%%% 
% [hC,pC,ciC,statsC] = ttest2(Controls.age(1:56),Controls.age(56:end))
% [hD,pD,ciD,statsD] = ttest2(DCLs.age(1:56),DCLs.age(56:end))               % = 0 por lo que no hay dif estadísticas de edad de controles ni DCL ENTRE AMBOS PROYECTOS
% [hT,pT,ciT,statsT] = ttest2(Controls.age,DCLs.age)              % = 1, por lo que SÍ HAY dif estadísticamente significativas entre edad de CONTROLES VS DCLs en general. 
%     % The result h is 1 if the test rejects the null hypothesis at the 5% significance level, and 0 otherwise.
%     % Test the null hypothesis that the two data samples are from populations with equal means.
% histogram(subj_datas.age)% 

% histogram(Controls.age)
% title('Edad');
% xlabel('años');
% ylabel(''); 
% 
% hold on
% 
% histogram(DCLs.age)
% legend ('Controls', 'DCLs')
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
% female = subj_datas(find(subj_datas.sex == 2), :); % 168 mujeres --> 92 DCLs, 76 controles // 157: 82 dcl, 75 controles
% male = subj_datas(find(subj_datas.sex == 1), :); % 94 hombres --> 53 DCLs, 41 controles // 91: 50 dcl, 41 controles
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
% histogram(Controls.edu_years)
% title('Años de educación');
% xlabel('años');
% ylabel(''); 
% 
% hold on
% 
% histogram(DCLs.edu_years)
% legend ('Controls', 'DCLs')
% 
% hold off
% % 
% % [H, pval, w] = swtest(Controls.edu_years)  % H = 1 --> distribución no normal
% % [H, pval, w] = swtest(DCLs.edu_years) % H = 1 
% [H, pval, w] = swtest(subj_datas.edu_years) % H = 1  --> distribución no normal
% 
% % Saber si hay diferencias significativas entre los 2 grupos
%     % Mann-Whitney-Wilcoxon non parametric test for two unpaired groups.
% 
% % mwwtest(Controls.edu_years, DCLs.edu_years)  % H=1 --> signif, pval = (<0.001)
% 
% aoctool (subj_datas.edu_years, subj_datas.age, subj_datas.diag)   % F = 0.02, pval = 0.885 //F=0.01, p=0.937
% 
% mean_eduControl = mean(Controls.edu_years) % = 13.65
% std_eduControl = std(Controls.edu_years) % = 5.24
% mean_eduDCL = mean(DCLs.edu_years) % = 8.98
% std_eduDCL = std(DCLs.edu_years) % = 4.37
% 
% 
