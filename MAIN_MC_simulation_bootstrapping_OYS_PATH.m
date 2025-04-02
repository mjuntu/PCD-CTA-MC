%% Compute MC simulation. Note: this version is with bootstrapping using the diagnostic path at the author's site. The other version estimates probabilities of misclassification from the patient data collection
% Mikael Brix, 14.3.2025
% Some of the parameters based on the following paper: 10.1016/j.jcct.2024.10.011
% Note in OYS, the treatment is directly to ICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define paths & user-defined variables
fileLoc = '/Users/mikaeljuntunen/Documents_Mikael/Work/Publications/Post-doc/PCD-CT-FINANCIALS/Codes_Git/Reclassification_data.xlsx';
N_CTA_patients = 1000;%465;%465; % Number of CTA patients a year, For OYS: 465 based on old statistics, 1000 based on actual need
cost_model = 'OYS'; % 'OYS' for our numbers, 'US' for reference publication. 
price_difference = 1.5e6; % Assume a 1.5 million euro difference between the price of Force and Alpha
switch cost_model
    case 'US'
        patient_limit = 15000; % In vecsey nagy
    case 'OYS'
        patient_limit = N_CTA_patients*10;
end
N_MC = 1e3; % Number of Monte Carlo iterations
p_level = 0.001; % Level of statistical significance
% For creating the result figures
fontsize = 14; 
fontname = 'Times new roman';
%% Define the critical fixed parameters that need to be defined
% CAD-RADS prevalences
%CAD_RADS_rates = [0.339,0.322,0.2,0.077,0.048,0.014]; % PROMISE TRIAL CAD-RADS 0, CAD-RADS 1, ...
CAD_RADS_rates = [0.487,0.303,0.116,0.073,0.021]; % CAD-RADS CAD-RADS 1, .... Note this modification was done to account for the fact that CAD-RADS 0 do not seek treatment

CAD_RADS_cumulative = cumsum(CAD_RADS_rates);
% Define functional test probabilities


prob_ICA_CAD_RADS4 = (0.618+0.648)/2;
%% Costs
switch cost_model
    case 'OYS' % costs in €
        % % 2024 numbers
        % costs.CTA = 1017; 
        % costs.ECG = 302; % stress ecg
        % costs.ECHO = 512 ; % Echocardiography
        % costs.SPECT = 1071; % Need to confirm
        % costs.ICA = 1501;
        % costs.ICA_PCI = 6149;
        % PCI_rate = 0.531; % 53.1 % of the time, a PCI is done for ICA patient
        % % y_lim = [0,5e6];
        % y_lim = [0,13e6]; % windowing limits
        % 2025 numbers
        costs.CTA = 1017; 
        costs.ECG = 302; % stress ecg
        costs.ECHO = 512 ; % Echocardiography
        costs.SPECT = 1071; % Need to confirm
        costs.ICA = 1591;
        costs.ICA_PCI = 6518;
        PCI_rate = 0.531; % 53.1 % of the time, a PCI is done for ICA patient
        y_lim = [0,15e6]; % windowing limits
    case 'US' % costs in $
        costs.CTA = 4819; 
        costs.ECG = 990; 
        costs.ECHO = 5310 ; 
        costs.SPECT = 3785; 
        costs.ICA = 83659;
        costs.ICA_PCI = 83659; % NOTE: They did not mention separate price for ICA vs ICA+PCI
        y_lim = [0,160e6];

end
%% Import the data
data = importfile(fileLoc,'Data');
%% Start the Monte Carlo for Force
CAD_RADS_SCORE_PREDS = zeros(patient_limit,N_MC);
CAD_RADS_SCORE_TRUES = zeros(patient_limit,N_MC);
treatment_costs = zeros(patient_limit,N_MC);

% Initialize total cost accumulators (outside parfor)
total_costs_temp = zeros(1, 4);  % Columns: ECG, ECHO, SPECT, ICA, ICA_PCI
test_counts_temp = zeros(1, 4);  

for ii_MC = 1:N_MC
    fprintf('Force MC round %d of %d\n',ii_MC,N_MC);

    % Define temporary accumulators for this Monte Carlo run
    total_costs_mc = zeros(1, 4); 
    test_counts_mc = zeros(1, 4); 

    treatment_costs_mc = zeros(patient_limit, 1);
    CAD_RADS_SCORE_PREDS_mc = zeros(patient_limit, 1);
    CAD_RADS_SCORE_TRUES_mc = zeros(patient_limit, 1);

    parfor ii_patient = 1:patient_limit
        % Sample a patient
        CAD_RADS_SCORE_TRUE = randsample(1:5, 1, true, CAD_RADS_rates);
        CAD_RADS_SCORE_TRUES_mc(ii_patient) = CAD_RADS_SCORE_TRUE;

        % Initialize local cost accumulators for this patient
        total_cost_local = zeros(1, 4); 
        test_count_local = zeros(1, 4); 

        % Find indices with the same true score
        idxs = find(data.QCACADRADS == CAD_RADS_SCORE_TRUE);
        if isempty(idxs)
            CAD_RADS_SCORE_PRED = 0;
            treatment_cost = costs.CTA;
        else
            % Perform bootstrapping from real data
            CAD_RADS_SCORE_PRED = data.ForceCADRADS(idxs(randsample(1:numel(idxs), 1)));
            treatment_cost = 0;
            ICA_cost = 0;

            switch CAD_RADS_SCORE_PRED
                case {0, 1, 2}
                    treatment_cost = costs.CTA;

                case {3,4,5}
                    if CAD_RADS_SCORE_TRUE<4
                        ICA_cost = costs.ICA;
                    else
                        ICA_cost = costs.ICA_PCI;
                    end
                    total_cost_local(4) = ICA_cost; % Index 5 = ICA_PCI
                    test_count_local(4) = 1;
                    treatment_cost = costs.CTA + ICA_cost;
            end
        end

        treatment_costs_mc(ii_patient) = treatment_cost;
        CAD_RADS_SCORE_PREDS_mc(ii_patient) = CAD_RADS_SCORE_PRED;

        % Use reduction variables to accumulate inside parfor
        total_costs_mc = total_costs_mc + total_cost_local;
        test_counts_mc = test_counts_mc + test_count_local;
    end

    % Accumulate results outside parfor
    treatment_costs(:,ii_MC) = treatment_costs_mc;
    CAD_RADS_SCORE_PREDS(:,ii_MC) = CAD_RADS_SCORE_PREDS_mc;
    CAD_RADS_SCORE_TRUES(:,ii_MC) = CAD_RADS_SCORE_TRUES_mc;
    
    total_costs_temp(ii_MC,:) = total_costs_mc;
    test_counts_temp(ii_MC,:) = test_counts_mc;
end

% Compute mean costs per test type
mean_costs_per_test = mean(total_costs_temp ./ (patient_limit),1); % Avoid division by zero
downstream_test_prevalence_Force = mean(test_counts_temp,1);
downstream_test_SD_Force = std(test_counts_temp,1);
% Display results
% Format the values with mean ± SD
formatted_values = arrayfun(@(mean_val, sd_val) sprintf('%.3f ± %.3f', mean_val, sd_val), ...
                            downstream_test_prevalence_Force, downstream_test_SD_Force, ...
                            'UniformOutput', false);

% Display the results in a formatted table
disp('Number of downstream tests for Force (Mean ± SD):');
disp(table({'ECG', 'ECHO', 'SPECT', 'ICA'}', formatted_values', ...
    'VariableNames', {'StudyType', 'Number of Examinations (Mean ± SD)'}));
cumulative_treatment_costs_Force = cumsum(treatment_costs,1);
cumulative_downstream_costs_Force = cumsum(treatment_costs-costs.CTA,1);

%% Start the Monte Carlo for Alpha
scores = [];
CAD_RADS_SCORES = zeros(patient_limit,N_MC);
CAD_RADS_SCORE_PREDS = zeros(patient_limit,N_MC);
CAD_RADS_SCORE_TRUES = zeros(patient_limit,N_MC);
treatment_costs = zeros(patient_limit,N_MC);

% Initialize total cost accumulators (outside parfor)
total_costs_temp = zeros(1, 4);  % Columns: ECG, ECHO, SPECT, ICA, ICA_PCI
test_counts_temp = zeros(1, 4);  

for ii_MC = 1:N_MC
    fprintf('Alpha MC round %d of %d\n',ii_MC,N_MC);

    % Define temporary accumulators for this Monte Carlo run
    total_costs_mc = zeros(1, 4); 
    test_counts_mc = zeros(1, 4); 

    treatment_costs_mc = zeros(patient_limit, 1);
    CAD_RADS_SCORE_PREDS_mc = zeros(patient_limit, 1);
    CAD_RADS_SCORE_TRUES_mc = zeros(patient_limit, 1);

    parfor ii_patient = 1:patient_limit
        % Sample a patient
        CAD_RADS_SCORE_TRUE = randsample(1:5, 1, true, CAD_RADS_rates);
        CAD_RADS_SCORE_TRUES_mc(ii_patient) = CAD_RADS_SCORE_TRUE;

        % Initialize local cost accumulators for this patient
        total_cost_local = zeros(1, 4); 
        test_count_local = zeros(1, 4); 

        % Find indices with the same true score
        idxs = find(data.QCACADRADS == CAD_RADS_SCORE_TRUE);
        if isempty(idxs)
            CAD_RADS_SCORE_PRED = 0;
            treatment_cost = costs.CTA;
        else
            % Perform bootstrapping from real data
            CAD_RADS_SCORE_PRED = data.AlphaCADRADS(idxs(randsample(1:numel(idxs), 1)));
            treatment_cost = 0;
            ICA_cost = 0;

            switch CAD_RADS_SCORE_PRED
                case {0, 1, 2}
                    treatment_cost = costs.CTA;

                case {3,4,5}
                    if CAD_RADS_SCORE_TRUE<4
                        ICA_cost = costs.ICA;
                    else
                        ICA_cost = costs.ICA_PCI;
                    end
                    total_cost_local(4) = ICA_cost; % Index 5 = ICA_PCI
                    test_count_local(4) = 1;
                    treatment_cost = costs.CTA + ICA_cost;
            end
        end

        treatment_costs_mc(ii_patient) = treatment_cost;
        CAD_RADS_SCORE_PREDS_mc(ii_patient) = CAD_RADS_SCORE_PRED;

        % Use reduction variables to accumulate inside parfor
        total_costs_mc = total_costs_mc + total_cost_local;
        test_counts_mc = test_counts_mc + test_count_local;
    end

    % Accumulate results **outside** parfor
    treatment_costs(:,ii_MC) = treatment_costs_mc;
    CAD_RADS_SCORE_PREDS(:,ii_MC) = CAD_RADS_SCORE_PREDS_mc;
    CAD_RADS_SCORE_TRUES(:,ii_MC) = CAD_RADS_SCORE_TRUES_mc;
    
    total_costs_temp(ii_MC,:) = total_costs_mc;
    test_counts_temp(ii_MC,:) = test_counts_mc;
end

% Compute mean costs per test type
mean_costs_per_test = mean(total_costs_temp ./ (patient_limit),1); % Avoid division by zero
downstream_test_prevalence_Alpha = mean(test_counts_temp,1);
downstream_test_SD_Alpha = std(test_counts_temp,1);
% Display results
% Format the values with mean ± SD
formatted_values = arrayfun(@(mean_val, sd_val) sprintf('%.3f ± %.3f', mean_val, sd_val), ...
                            downstream_test_prevalence_Alpha, downstream_test_SD_Alpha, ...
                            'UniformOutput', false);

% Display the results in a formatted table
disp('Number of downstream tests for Alpha (Mean ± SD):');
disp(table({'ECG', 'ECHO', 'SPECT', 'ICA'}', formatted_values', ...
    'VariableNames', {'StudyType', 'Number of Examinations (Mean ± SD)'}));
cumulative_treatment_costs_Alpha = cumsum(treatment_costs,1);
cumulative_downstream_costs_Alpha = cumsum(treatment_costs-costs.CTA,1);
%% Visualization of downstream costs considering the price difference
cumulative_CT_cost = (1:size(cumulative_downstream_costs_Force,1))*price_difference/size(cumulative_downstream_costs_Force,1);

% Assume CAD_RADS_SCORE_PREDS is an (N_patients x N_MC_iterations) matrix
mean_preds = mean(cumulative_downstream_costs_Force, 2); % Mean over MC iterations
std_preds = std(cumulative_downstream_costs_Force, 0, 2); % Standard deviation over MC iterations

ci_upper = mean_preds + 1.96 * std_preds; % Upper 95% CI
ci_lower = mean_preds - 1.96 * std_preds; % Lower 95% CI

% Define x-axis (patients)
x = (1:size(cumulative_downstream_costs_Force, 1))';

% Plot the mean prediction
figure;
plot(x, mean_preds, 'r', 'LineWidth', 2); % Mean line

hold on;
% Fill area for confidence interval
fill([x; flipud(x)], [ci_upper; flipud(ci_lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Labels and title
xlabel('Number of patients');
ylabel('Cumulative downstream costs (with price difference) (€)');
%title('Cumulative downstream costs accounting for the higher price of PCD-CT');
set(gca,'fontsize',fontsize,'fontname',fontname);
% Visualization Alpha
% Assume CAD_RADS_SCORE_PREDS is an (N_patients x N_MC_iterations) matrix
%mean_preds = mean(cumulative_downstream_costs_Alpha, 2)+cumulative_CT_cost'; % Mean over MC iterations
mean_preds = mean(cumulative_downstream_costs_Alpha, 2)+cumulative_CT_cost(end)'; % Mean over MC iterations

std_preds = std(cumulative_downstream_costs_Alpha, 0, 2); % Standard deviation over MC iterations


ci_upper = mean_preds + 1.96 * std_preds; % Upper 95% CI
ci_lower = mean_preds - 1.96 * std_preds; % Lower 95% CI
% Determine the significant number of patients for cost difference
p_vals = zeros(size(cumulative_treatment_costs_Alpha,1),1);
for ii = 1:size(cumulative_treatment_costs_Alpha,1)
   % [h,p] = ttest2(cumulative_treatment_costs_Alpha(ii,:)+cumulative_CT_cost(end),cumulative_treatment_costs_Force(ii,:),'Alpha',p_level);
   [p, h] = ranksum(cumulative_treatment_costs_Alpha(ii,:)+cumulative_CT_cost(end), cumulative_treatment_costs_Force(ii,:)); 
   % figure;boxplot([cumulative_treatment_costs_Alpha(ii,:)'+cumulative_CT_cost(end),cumulative_treatment_costs_Force(ii,:)']) 
   p_vals(ii,1) = p;
    if p<p_level && mean(cumulative_treatment_costs_Alpha(ii,:)+cumulative_CT_cost(end)) < mean(cumulative_treatment_costs_Force(ii,:))
        fprintf('\nStatistically significant difference reached at iteration %d\n',ii);
        break;
    end
end
% Define x-axis (patients)
x = (1:size(cumulative_downstream_costs_Alpha, 1))';

% Plot the mean prediction
plot(x, mean_preds, 'b', 'LineWidth', 2); % Mean line

hold on;
% Fill area for confidence interval
fill([x; flipud(x)], [ci_upper; flipud(ci_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Labels and title
grid on;
hold off;
legend('EID-CT','','PCD-CT')
ylim(y_lim)


%% Visualization of cumulative downstream costs (no price difference) 

% Assume CAD_RADS_SCORE_PREDS is an (N_patients x N_MC_iterations) matrix
mean_preds = mean(cumulative_downstream_costs_Force, 2); % Mean over MC iterations
std_preds = std(cumulative_downstream_costs_Force, 0, 2); % Standard deviation over MC iterations
mean_pred_Force = mean_preds(end);
std_pred_Force = std_preds(end);

ci_upper = mean_preds + 1.96 * std_preds; % Upper 95% CI
ci_lower = mean_preds - 1.96 * std_preds; % Lower 95% CI

% Define x-axis (patients)
x = (1:size(cumulative_downstream_costs_Force, 1))';

% Plot the mean prediction
figure;
plot(x, mean_preds, 'r', 'LineWidth', 2); % Mean line

hold on;
% Fill area for confidence interval
fill([x; flipud(x)], [ci_upper; flipud(ci_lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Labels and title
xlabel('Number of patients');
ylabel('Cumulative treatment cost (€)');
% title('Cumulative downstream costs');
set(gca,'fontsize',fontsize,'fontname',fontname);

% Visualization Alpha
% Assume CAD_RADS_SCORE_PREDS is an (N_patients x N_MC_iterations) matrix
mean_preds = mean(cumulative_downstream_costs_Alpha, 2); % Mean over MC iterations
std_preds = std(cumulative_downstream_costs_Alpha, 0, 2); % Standard deviation over MC iterations
mean_pred_Alpha = mean_preds(end);
std_pred_Alpha = std_preds(end);

ci_upper = mean_preds + 1.96 * std_preds; % Upper 95% CI
ci_lower = mean_preds - 1.96 * std_preds; % Lower 95% CI

% Define x-axis (patients)
x = (1:size(cumulative_downstream_costs_Alpha, 1))';

% Plot the mean prediction
plot(x, mean_preds, 'b', 'LineWidth', 2); % Mean line

hold on;
% Fill area for confidence interval
fill([x; flipud(x)], [ci_upper; flipud(ci_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Labels and title
grid on;
hold off;
legend('EID-CT','','PCD-CT')
ylim(y_lim)
%% Mean downstream costs
fprintf('\nMean downstream costs:\n')
fprintf('Force: %4.3f ± %4.3f\n',mean_pred_Force/patient_limit,std_pred_Force/patient_limit)
fprintf('Alpha: %4.3f ± %4.3f\n',mean_pred_Alpha/patient_limit,std_pred_Alpha/patient_limit)
fprintf('Difference: %4.3f\n',mean_pred_Force/patient_limit-mean_pred_Alpha/patient_limit)
% fprintf('\nReduction in functional tests %4.2f\n',(functional_tests_Force-functional_tests_Alpha)/functional_tests_Force*100)
% fprintf('\nReduction in ICAs %4.2f\n',(ICAs_Alpha-ICAs_Force)/ICAs_Force*100)


%% Visualization of total costs

% Assume CAD_RADS_SCORE_PREDS is an (N_patients x N_MC_iterations) matrix
mean_preds = mean(cumulative_treatment_costs_Force, 2); % Mean over MC iterations
std_preds = std(cumulative_treatment_costs_Force, 0, 2); % Standard deviation over MC iterations
mean_pred_Force = mean_preds(end);
std_pred_Force = std_preds(end);

ci_upper = mean_preds + 1.96 * std_preds; % Upper 95% CI
ci_lower = mean_preds - 1.96 * std_preds; % Lower 95% CI

% Define x-axis (patients)
x = (1:size(cumulative_treatment_costs_Force, 1))';

% Plot the mean prediction
figure;
plot(x, mean_preds, 'r', 'LineWidth', 2); % Mean line

hold on;
% Fill area for confidence interval
fill([x; flipud(x)], [ci_upper; flipud(ci_lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Labels and title
xlabel('Number of patients');
ylabel('Cumulative total cost (€)');
title('Cumulative total costs');

% Visualization Alpha
% Assume CAD_RADS_SCORE_PREDS is an (N_patients x N_MC_iterations) matrix
mean_preds = mean(cumulative_treatment_costs_Alpha, 2); % Mean over MC iterations
std_preds = std(cumulative_treatment_costs_Alpha, 0, 2); % Standard deviation over MC iterations
mean_pred_Alpha = mean_preds(end);
std_pred_Alpha = std_preds(end);

ci_upper = mean_preds + 1.96 * std_preds; % Upper 95% CI
ci_lower = mean_preds - 1.96 * std_preds; % Lower 95% CI

% Define x-axis (patients)
x = (1:size(cumulative_treatment_costs_Alpha, 1))';

% Plot the mean prediction
plot(x, mean_preds, 'b', 'LineWidth', 2); % Mean line

hold on;
% Fill area for confidence interval
fill([x; flipud(x)], [ci_upper; flipud(ci_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Labels and title
grid on;
hold off;
legend('Force','','Alpha')
ylim(y_lim)
%% Mean total costs
fprintf('\nMean total costs:\n')
fprintf('Force: %4.3f ± %4.3f\n',mean_pred_Force/patient_limit,std_pred_Force/patient_limit)
fprintf('Alpha: %4.3f ± %4.3f\n',mean_pred_Alpha/patient_limit,std_pred_Alpha/patient_limit)
fprintf('Difference: %4.3f\n',mean_pred_Force/patient_limit-mean_pred_Alpha/patient_limit)
%fprintf('\nReduction in functional tests %4.2f\n',(functional_tests_Force-functional_tests_Alpha)/functional_tests_Force*100)
%fprintf('\nReduction in ICAs %4.2f\n',(ICAs_Alpha-ICAs_Force)/ICAs_Force*100)
%% Determine the significant number of patients for cost difference
%% Determine the significant number of patients for cost difference
p_vals = zeros(size(cumulative_treatment_costs_Alpha,1),1);
for ii = 1:size(cumulative_treatment_costs_Alpha,1)
    [h,p] = ttest(cumulative_treatment_costs_Alpha(ii,:),cumulative_treatment_costs_Force(ii,:));
    p_vals(ii,1) = p;
    % if p<p_level
    %     fprintf('\nStatistically significant difference reached at iteration %d\n',ii);
    %     break;
    % end
end
%%
figure;semilogy(1:numel(p_vals),p_vals);ylim([0,0.1]);

%% Number of downstream examinations
% Display results
% Display results
% Format the values with mean ± SD
formatted_values = arrayfun(@(mean_val, sd_val) sprintf('%.3f ± %.3f', mean_val, sd_val), ...
                            downstream_test_prevalence_Force, downstream_test_SD_Force, ...
                            'UniformOutput', false);

% Display the results in a formatted table
disp('Number of downstream tests for Force (Mean ± SD):');
disp(table({'ECG', 'ECHO', 'SPECT', 'ICA'}', formatted_values', ...
    'VariableNames', {'StudyType', 'Number of Examinations (Mean ± SD)'}));

formatted_values = arrayfun(@(mean_val, sd_val) sprintf('%.3f ± %.3f', mean_val, sd_val), ...
                            downstream_test_prevalence_Alpha, downstream_test_SD_Alpha, ...
                            'UniformOutput', false);

% Display the results in a formatted table
disp('Number of downstream tests for Alpha (Mean ± SD):');
disp(table({'ECG', 'ECHO', 'SPECT', 'ICA'}', formatted_values', ...
    'VariableNames', {'StudyType', 'Number of Examinations (Mean ± SD)'}));
