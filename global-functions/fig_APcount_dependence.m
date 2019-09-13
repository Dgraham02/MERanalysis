% Global Function 
% Purpose: Compare AP aligned LFP RMS to AP count 
% Notes: 

function [f] = fig_APcount_dependence(general_data, post_hp_data)

patient_path = general_data.patient_path;
patient_id = general_data.patient_id;

fontsize = 14;

f = figure('Name', 'AP Count vs RMS');clf;
% Plot RMS Values 
subplot(1,3,1)
for i = 1:length(blocks)
    plot(flipud(rmsMat(:,i)),'o-')
    hold on
end
title('LFP Average RMS','fontsize',fontsize/1.2)
legend('140Hz','20Hz','250Hz','70Hz','Location','northwest')
ylabel('Voltage (uV)','fontsize',fontsize/1.2)
xticks([1 2 3])
xticklabels({'Baseline','OFF','ON'})

% Plot AP Count Values 
subplot(1,3,2)
for i = 1:length(blocks)
    plot(flipud(countMat(:,i)),'o-')
    hold on
end
title('AP Count','fontsize',fontsize/1.2)
legend('140Hz','20Hz','250Hz','70Hz','Location','northeast')
ylabel('# of Action Potentials','fontsize',fontsize/1.2)
xticks([1 2 3])
xticklabels({'Baseline','OFF','ON'})

% Plot AP Count vs RMS 
xcount = reshape(countMat,[],1);
xrms = reshape(rmsMat,[],1);

subplot(1,3,3)
for i = 1:length(blocks)
    scatter(flipud(countMat(:,i)),flipud(rmsMat(:,i)))
    hold on
end
%scatter(xcount,xrms)
title('AP Count vs RMS','fontsize',fontsize/1.2)
sgtitle('AP Count & RMS Dependence: Contralateral Patient #1','fontsize',fontsize*1.2)
xlabel('AP Count','fontsize',fontsize/1.2)
ylabel('Voltage (uV)','fontsize',fontsize/1.2)
legend('140Hz','20Hz','250Hz','70Hz','Location','northeast')