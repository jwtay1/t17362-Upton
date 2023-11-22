%Analyze data
clearvars
clc

exp = load('experiment.mat');
control = load('control.mat');

%Calculate the change in DAPI intensity?
expDAPInorm = zeros(numel(exp.combinedCellData), 3);
for iCell = 1:numel(exp.combinedCellData)

    expDAPInorm(iCell, :) = exp.combinedCellData(iCell).DAPIbleach / (exp.combinedCellData(iCell).DAPIbleach(1));


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
DAPIchange = [ctlDAPInorm(:, 3); expDAPInorm(:, 3)];
grouping = [ones(size(ctlDAPInorm, 1), 1); ...
    ones(size(expDAPInorm, 1), 1) * 2];

[P, t, st] = anova1(DAPIchange, grouping);
%multcompare(st)

