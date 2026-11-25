%% Function for extracting p-values of pairwise comparisons from LME output

% Arguments: 
% lme: LinearMixedModel object from fitlme
% groupVar: name of grouping variable (string), e.g., 'Group'

function pairwisePvals = lme_pairwise_pvals(lme, groupVar)

    tbl = lme.Variables;
    groups = categories(categorical(tbl.(groupVar)));
    coefNames = lme.CoefficientNames;

    results = [];
    row = 1;

    for i = 1:length(groups)
        for j = i+1:length(groups)

            % Build contrast vector:
            %   +1 for Group_i coefficient
            %   -1 for Group_j coefficient
            contrast = zeros(1, length(coefNames));

            coef_i = [groupVar '_' groups{i}];
            coef_j = [groupVar '_' groups{j}];

            if any(strcmp(coefNames, coef_i))
                contrast(strcmp(coefNames, coef_i)) = 1;
            end

            if any(strcmp(coefNames, coef_j))
                contrast(strcmp(coefNames, coef_j)) = -1;
            end

            [p, F, df1, df2] = coefTest(lme, contrast);
            
            % Compute estimate and t-statistic 
            beta = lme.Coefficients.Estimate';
            estimate = dot(contrast, beta);
            se = sqrt(contrast * lme.CoefficientCovariance * contrast');
            tstat = estimate / se;

            % Store row
            results{row,1} = groups{i};
            results{row,2} = groups{j};
            results{row,3} = estimate;
            results{row,4} = tstat;
            results{row,5} = p;
            results{row,6} = df1;
            results{row,7} = df2;

            row = row + 1;
        end
    end

    % Convert to table
    pairwisePvals = cell2table(results, ...
        'VariableNames', {'Group1','Group2','Estimate','t','p','DF1','DF2'});
end