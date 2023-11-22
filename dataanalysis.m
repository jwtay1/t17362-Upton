%Analyze data
clearvars
clc

exp = load('experiment.mat');
exp2 = load('experiment2.mat');
control = load('control.mat');

%Calculate the change in DAPI intensity?
expDAPInorm = zeros(numel(exp.combinedCellData), 3);
for iCell = 1:numel(exp.combinedCellData)

    expDAPInorm(iCell, :) = exp.combinedCellData(iCell).DAPIbleach / (exp.combinedCellData(iCell).DAPIbleach(1));


end

%Calculate the change in DAPI intensity?
expDAPInorm2 = zeros(numel(exp2.combinedCellData), 3);
for iCell = 1:numel(exp2.combinedCellData)

    expDAPInorm2(iCell, :) = exp2.combinedCellData(iCell).DAPIbleach / (exp2.combinedCellData(iCell).DAPIbleach(1));


end


%Calculate the change in DAPI intensity?
ctlDAPInorm = zeros(numel(control.combinedCellData), 3);
for iCell = 1:numel(control.combinedCellData)

    ctlDAPInorm(iCell, :) = control.combinedCellData(iCell).DAPIbleach / (control.combinedCellData(iCell).DAPIbleach(1));


end

histogram(expDAPInorm(:, 3))
hold on
histogram(ctlDAPInorm(:, 3))
hold off
%%
%Compile data into two grpuds
DAPIchange = [ctlDAPInorm(:, 3); expDAPInorm(:, 3); expDAPInorm2(:, 3)];
grouping = [ones(size(ctlDAPInorm, 1), 1); ...
    ones(size(expDAPInorm, 1), 1) * 2; ...
    ones(size(expDAPInorm2, 1), 1) * 3];

[P, t, st] = anova1(DAPIchange, grouping);
%multcompare(st)

