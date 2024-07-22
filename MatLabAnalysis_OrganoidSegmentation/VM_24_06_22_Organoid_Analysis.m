clc
close all
clear all
%%

FolderName=% add folder path with the excel files from OrgaQuant
folder=strcat('');
files = dir(folder);

data=struct([]);
maxVols=zeros(1,36);

for j=1:36
    k=3*j+1;
    name=files(k).name;
    FileName=strcat(folder,'\',name);
    name=name(1:end-9);
    data(j).name=name;
    data(j).case=str2double(name(1));
    data(j).rep=str2double(name(end));
    data(j).day=str2double(name(6:7));
    
    matrix=[];
    matrix=table2array(readtable(FileName,'NumHeaderLines',1));
    matrix=matrix(:,5:6);

    matrix1=[];
    matrix1(:,1)=(matrix(:,1)+matrix(:,2))/2;
    matrix1(matrix1>800)=[];

    matrixArea=[];
    matrixArea(:,1)=round(pi*(matrix1(:,1)/2).^2);
    
    maxVols(j)=max(matrixArea(:));
    data(j).areaabs=matrixArea;
    data(j).arearel=matrixArea/sum(matrixArea(:));
    data(j).areatot=sum(matrixArea(:))/(1080*1920);

    [ycdf,xcdf]=ecdf(matrixArea);
    data(j).cdf=[xcdf ycdf];

    %[ycdfrel,xcdfrel]=ecdf(matrixArea/sum(matrixArea(:)));
    [ycdfrel,xcdfrel]=ecdf(matrixArea/max(matrixArea(:)));
    data(j).cdfrel=[xcdfrel,ycdfrel];
    
    maxsize=408849;
    % maxsize=89462;
    compL=200-size(xcdf,1);
    xend=round(xcdf(end));
    compx=linspace(xend,maxsize,compL)';
    compy=ones(compL,1);
    
    cdfc=[xcdf ycdf;
          compx compy];

    data(j).cdfc=cdfc;
    data(j).cdfcnorm=[cdfc(:,1)/maxsize cdfc(:,2)];  
    data(j).AUC=trapz(cdfc(:,1)/maxsize,cdfc(:,2));
    data(j).Inv_AUC=1-data(j).AUC;
end


%% Statistical tests

Medium=[1 1 1 2 2 2 3 3 3]';

A_tot_7=[data(16:18).areatot data(25:27).areatot data(34:36).areatot]';
Inv_AUC_7=[data(16:18).Inv_AUC data(25:27).Inv_AUC data(34:36).Inv_AUC]';

A_tot_12=[data(10:12).areatot data(19:21).areatot data(28:30).areatot]';
Inv_AUC_12=[data(10:12).Inv_AUC data(19:21).Inv_AUC data(28:30).Inv_AUC]';

A_tot_19=[data(13:15).areatot data(22:24).areatot data(31:33).areatot]';
Inv_AUC_19=[data(13:15).Inv_AUC data(22:24).Inv_AUC data(31:33).Inv_AUC]';

[p_A7, tbl_A7, stats_A7] = anova1(A_tot_7, Medium,'off');
c_A7 = multcompare(stats_A7, 'CType', 'tukey-kramer','Display', 'off');
[p_AUC7, tbl_AUC7, stats_AUC7] = anova1(Inv_AUC_7, Medium,'off');
c_AUC7 = multcompare(stats_AUC7, 'CType', 'tukey-kramer','Display', 'off');

[p_A12, tbl_A12, stats_A12] = anova1(A_tot_12, Medium,'off');
c_A12 = multcompare(stats_A12, 'CType', 'tukey-kramer','Display', 'off');
[p_AUC12, tbl_AUC12, stats_AUC12] = anova1(Inv_AUC_12, Medium,'off');
c_AUC12 = multcompare(stats_AUC12, 'CType', 'tukey-kramer','Display', 'off');

[p_A19, tbl_A19, stats_A19] = anova1(A_tot_19, Medium,'off');
c_A19 = multcompare(stats_A19, 'CType', 'tukey-kramer','Display', 'off');
[p_AUC19, tbl_AUC19, stats_AUC19] = anova1(Inv_AUC_19, Medium,'off');
c_AUC19 = multcompare(stats_AUC19, 'CType', 'tukey-kramer','Display', 'off');
%% Plot total area

Areatot=zeros(1,36);
for k=1:36
    Areatot(1,k)=data(k).areatot;
end

meanArea=100*[mean(Areatot(16:18)) mean(Areatot(25:27)) mean(Areatot(34:36));
              mean(Areatot(10:12)) mean(Areatot(19:21)) mean(Areatot(28:30))
              mean(Areatot(13:15)) mean(Areatot(22:24)) mean(Areatot(31:33))];

stdArea=(200/sqrt(3))*[std(Areatot(16:18)) std(Areatot(25:27)) std(Areatot(34:36));
                       std(Areatot(10:12)) std(Areatot(19:21)) std(Areatot(28:30))
                       std(Areatot(13:15)) std(Areatot(22:24)) std(Areatot(31:33))];

close all
figure(1);
hold on
box on
x = [1,2,3];
b = bar(x,meanArea,'LineWidth',1.5);
xCnt=cell2mat(get(b,'XEndPoints'))'; 
errorbar(xCnt, (meanArea),stdArea, 'k', 'LineStyle','none','LineWidth',1.5);
yCnt = meanArea + stdArea; 
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18,'FontWeight','Bold',  'LineWidth', 2);
ylabel('Tracheosphere Area/Image Area (%)');
ylim([0, 80]);
set(gca,'xtick',1:3,'xticklabel',{'Day 7','Day 12','Day 19'})
b(1).FaceColor = [0.2 0.2 0.2];
b(2).FaceColor = [0.4 0.4 0.4];
b(3).FaceColor = [0.6 0.6 0.6];
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18,'FontWeight','Bold',  'LineWidth', 2);
% Adding text annotations on specified bars
% Add statistical comparison bars
bar_offset = 2;  % Distance above the error bars
line([xCnt(2,1), xCnt(2,2)], [yCnt(2,3) + bar_offset, yCnt(2,3) + bar_offset], 'Color', 'k', 'LineWidth', 1.5)
line([xCnt(2,1), xCnt(2,3)], [yCnt(2,3) + bar_offset+4, yCnt(2,3) + bar_offset+4], 'Color', 'k', 'LineWidth', 1.5)
line([xCnt(2,2), xCnt(2,3)], [yCnt(2,3) + bar_offset+2, yCnt(2,3) + bar_offset+2], 'Color', 'k', 'LineWidth', 1.5)
text((xCnt(2,1)+xCnt(2,2))/2, yCnt(2,3) + bar_offset + 0.5, '***', 'HorizontalAlignment', 'center', 'FontSize', 18)
text((xCnt(2,1)+xCnt(2,3))/2, yCnt(2,3) + 4+ bar_offset + 0.5, '***', 'HorizontalAlignment', 'center', 'FontSize', 18)
text((xCnt(2,2)+xCnt(2,3))/2, yCnt(2,3) + 2+ bar_offset + 0.5, '*', 'HorizontalAlignment', 'center', 'FontSize', 18)

%line([xCnt(3,1), xCnt(3,2)], [yCnt(3,3) + bar_offset, yCnt(3,3) + bar_offset], 'Color', 'k', 'LineWidth', 1.5)
line([xCnt(3,1), xCnt(3,3)], [yCnt(3,3) + bar_offset+4, yCnt(3,3) + bar_offset+4], 'Color', 'k', 'LineWidth', 1.5)
line([xCnt(3,2), xCnt(3,3)], [yCnt(3,3) + bar_offset+2, yCnt(3,3) + bar_offset+2], 'Color', 'k', 'LineWidth', 1.5)
%text((xCnt(3,1)+xCnt(3,2))/2, yCnt(3,3) + bar_offset + 0.5, '***', 'HorizontalAlignment', 'center', 'FontSize', 18)
text((xCnt(3,1)+xCnt(3,3))/2, yCnt(3,3) + 4+ bar_offset + 0.5, '*', 'HorizontalAlignment', 'center', 'FontSize', 18)
text((xCnt(3,2)+xCnt(3,3))/2, yCnt(3,3) + 2+ bar_offset + 0.5, '*', 'HorizontalAlignment', 'center', 'FontSize', 18)

legend(b,'Medium 1','Medium 2','Medium 3', 'location','northwest');

%% CDF normalized largest organoid in the dataset

figure(2)

subplot(231)
box on
hold on

plot(data(16).cdfcnorm(:,1),data(16).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(25).cdfcnorm(:,1),data(25).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(34).cdfcnorm(:,1),data(34).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(data(17).cdfcnorm(:,1),data(17).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(18).cdfcnorm(:,1),data(18).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(26).cdfcnorm(:,1),data(26).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(27).cdfcnorm(:,1),data(27).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(35).cdfcnorm(:,1),data(35).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(data(36).cdfcnorm(:,1),data(36).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(0:0.005:0.5,0.8*ones(1,101),'--k','LineWidth',1.5);
plot(0.5*ones(1,41),0.8:0.005:1,'--k','LineWidth',1.5);
%legend('Media 1','Media 2','Media 3', 'location','southeast');
xlabel('Normalized Tracheosphere Area');
ylabel('CDF');
title('Day 7');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18,'FontWeight','Bold',  'LineWidth', 2);
hold off

subplot(232)
box on
hold on
plot(data(10).cdfcnorm(:,1),data(10).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(19).cdfcnorm(:,1),data(19).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(28).cdfcnorm(:,1),data(28).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(data(11).cdfcnorm(:,1),data(11).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(12).cdfcnorm(:,1),data(12).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(20).cdfcnorm(:,1),data(20).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(21).cdfcnorm(:,1),data(21).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(29).cdfcnorm(:,1),data(29).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(data(30).cdfcnorm(:,1),data(30).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(0:0.005:0.5,0.8*ones(1,101),'--k','LineWidth',1.5);
plot(0.5*ones(1,41),0.8:0.005:1,'--k','LineWidth',1.5);
%legend('Media 1','Media 2','Media 3', 'location','southeast');
xlabel('Normalized Tracheosphere Area');
ylabel('CDF');
title('Day 12');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18,'FontWeight','Bold',  'LineWidth', 2);
hold off

subplot(233)
box on
hold on
plot(data(13).cdfcnorm(:,1),data(13).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(22).cdfcnorm(:,1),data(22).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(31).cdfcnorm(:,1),data(31).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(data(14).cdfcnorm(:,1),data(14).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(15).cdfcnorm(:,1),data(15).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(23).cdfcnorm(:,1),data(23).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(24).cdfcnorm(:,1),data(24).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(32).cdfcnorm(:,1),data(32).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(data(33).cdfcnorm(:,1),data(33).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(0:0.005:0.5,0.8*ones(1,101),'--k','LineWidth',1.5);
plot(0.5*ones(1,41),0.8:0.005:1,'--k','LineWidth',1.5);
%legend('Media 1','Media 2','Media 3', 'location','southeast');
xlabel('Normalized Tracheosphere Area');
ylabel('CDF');
title('Day 19');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18,'FontWeight','Bold',  'LineWidth', 2);
hold off


subplot(234)
box on
hold on

plot(data(16).cdfcnorm(:,1),data(16).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(25).cdfcnorm(:,1),data(25).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(34).cdfcnorm(:,1),data(34).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(data(17).cdfcnorm(:,1),data(17).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(18).cdfcnorm(:,1),data(18).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(26).cdfcnorm(:,1),data(26).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(27).cdfcnorm(:,1),data(27).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(35).cdfcnorm(:,1),data(35).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(data(36).cdfcnorm(:,1),data(36).cdfcnorm(:,2),'-r','LineWidth',1.5);
%legend('Media 1','Media 2','Media 3', 'location','southeast');
xlabel('Normalized Tracheosphere Area');
ylabel('CDF');
title('Day 7');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18,'FontWeight','Bold',  'LineWidth', 2);
xlim([0 0.5]);
ylim([0.8 1]);
hold off

subplot(235)
box on
hold on
plot(data(10).cdfcnorm(:,1),data(10).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(19).cdfcnorm(:,1),data(19).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(28).cdfcnorm(:,1),data(28).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(data(11).cdfcnorm(:,1),data(11).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(12).cdfcnorm(:,1),data(12).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(20).cdfcnorm(:,1),data(20).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(21).cdfcnorm(:,1),data(21).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(29).cdfcnorm(:,1),data(29).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(data(30).cdfcnorm(:,1),data(30).cdfcnorm(:,2),'-r','LineWidth',1.5);
legend('Medium 1','Medium 2','Medium 3', 'location','southoutside','Orientation','Horizontal');
xlabel('Normalized Tracheosphere Area');
ylabel('CDF');
title('Day 12');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18,'FontWeight','Bold',  'LineWidth', 2);
xlim([0 0.5]);
ylim([0.8 1]);
hold off

subplot(236)
box on
hold on
plot(data(13).cdfcnorm(:,1),data(13).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(22).cdfcnorm(:,1),data(22).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(31).cdfcnorm(:,1),data(31).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(data(14).cdfcnorm(:,1),data(14).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(15).cdfcnorm(:,1),data(15).cdfcnorm(:,2),'-b','LineWidth',1.5);
plot(data(23).cdfcnorm(:,1),data(23).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(24).cdfcnorm(:,1),data(24).cdfcnorm(:,2),'-g','LineWidth',1.5);
plot(data(32).cdfcnorm(:,1),data(32).cdfcnorm(:,2),'-r','LineWidth',1.5);
plot(data(33).cdfcnorm(:,1),data(33).cdfcnorm(:,2),'-r','LineWidth',1.5);
xlabel('Normalized Tracheosphere Area');
ylabel('CDF');
title('Day 19');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18,'FontWeight','Bold',  'LineWidth', 2);
xlim([0 0.5]);
ylim([0.8 1]);
hold off

%% Plot AUC

AUCtot=zeros(1,36);
for k=1:36
    AUCtot(1,k)=data(k).AUC;
end

AUCtot=1-AUCtot;

meanAUC=[mean(AUCtot(16:18)) mean(AUCtot(25:27)) mean(AUCtot(34:36));
             mean(AUCtot(10:12)) mean(AUCtot(19:21)) mean(AUCtot(28:30))
             mean(AUCtot(13:15)) mean(AUCtot(22:24)) mean(AUCtot(31:33))];

AUCAll=zeros(3,27);
    

stdAUC=(2/sqrt(3))*[std(AUCtot(16:18)) std(AUCtot(25:27)) std(AUCtot(34:36));
                      std(AUCtot(10:12)) std(AUCtot(19:21)) std(AUCtot(28:30))
                      std(AUCtot(13:15)) std(AUCtot(22:24)) std(AUCtot(31:33))];

figure(3);
hold on
box on
x = [1,2,3];
b = bar(x,meanAUC,'LineWidth',1.5);
xCnt=cell2mat(get(b,'XEndPoints'))'; 
yCnt = (meanAUC) + stdAUC; 
errorbar(xCnt, (meanAUC),stdAUC, 'k', 'LineStyle','none','LineWidth',1.5);
ylabel('1-AUC');
ylim([0 0.31]);
set(gca,'xtick',1:3,'xticklabel',{'Day 7','Day 12','Day 19'})
b(1).FaceColor = [0.2 0.2 0.2];
b(2).FaceColor = [0.4 0.4 0.4];
b(3).FaceColor = [0.6 0.6 0.6];
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18,'FontWeight','Bold',  'LineWidth', 2);


bar_offset = 0.02;  % Distance above the error bars

%line([xCnt(1,1), xCnt(1,2)], [yCnt(1,3) + bar_offset, yCnt(1,3) + bar_offset], 'Color', 'k', 'LineWidth', 1.5)
line([xCnt(1,1), xCnt(1,3)], [yCnt(1,3) + bar_offset, yCnt(1,3) + bar_offset], 'Color', 'k', 'LineWidth', 1.5)
%line([xCnt(1,2), xCnt(1,3)], [yCnt(1,3) + bar_offset+0.02, yCnt(1,3) + bar_offset+0.02], 'Color', 'k', 'LineWidth', 1.5)
%text((xCnt(1,1)+xCnt(1,2))/2, yCnt(1,3) + bar_offset +0.005, '***', 'HorizontalAlignment', 'center', 'FontSize', 18)
text((xCnt(1,1)+xCnt(1,3))/2, yCnt(1,3) + bar_offset +0.005, '*', 'HorizontalAlignment', 'center', 'FontSize', 18)
%text((xCnt(1,2)+xCnt(1,3))/2, yCnt(1,3) + 2+ bar_offset +0.005, '*', 'HorizontalAlignment', 'center', 'FontSize', 18)

line([xCnt(2,1), xCnt(2,2)], [yCnt(2,3) + bar_offset, yCnt(2,3) + bar_offset], 'Color', 'k', 'LineWidth', 1.5)
line([xCnt(2,1), xCnt(2,3)], [yCnt(2,3) + bar_offset+0.04, yCnt(2,3) + bar_offset+0.04], 'Color', 'k', 'LineWidth', 1.5)
line([xCnt(2,2), xCnt(2,3)], [yCnt(2,3) + bar_offset+0.02, yCnt(2,3) + bar_offset+0.02], 'Color', 'k', 'LineWidth', 1.5)
text((xCnt(2,1)+xCnt(2,2))/2, yCnt(2,3) + bar_offset +0.005, '***', 'HorizontalAlignment', 'center', 'FontSize', 18)
text((xCnt(2,1)+xCnt(2,3))/2, yCnt(2,3) +0.04+ bar_offset +0.005, '***', 'HorizontalAlignment', 'center', 'FontSize', 18)
text((xCnt(2,2)+xCnt(2,3))/2, yCnt(2,3) +0.02+ bar_offset +0.005, '***', 'HorizontalAlignment', 'center', 'FontSize', 18)

%line([xCnt(3,1), xCnt(3,2)], [yCnt(3,3) + bar_offset, yCnt(3,3) + bar_offset], 'Color', 'k', 'LineWidth', 1.5)
line([xCnt(3,1), xCnt(3,3)], [yCnt(3,3) + bar_offset+0.04, yCnt(3,3) + bar_offset+0.04], 'Color', 'k', 'LineWidth', 1.5)
line([xCnt(3,2), xCnt(3,3)], [yCnt(3,3) + bar_offset+0.02, yCnt(3,3) + bar_offset+0.02], 'Color', 'k', 'LineWidth', 1.5)
%text((xCnt(3,1)+xCnt(3,2))/2, yCnt(3,3) + bar_offset +0.005, '***', 'HorizontalAlignment', 'center', 'FontSize', 18)
text((xCnt(3,1)+xCnt(3,3))/2, yCnt(3,3) +0.04+ bar_offset +0.005, '*', 'HorizontalAlignment', 'center', 'FontSize', 18)
text((xCnt(3,2)+xCnt(3,3))/2, yCnt(3,3) +0.02+ bar_offset +0.005, '*', 'HorizontalAlignment', 'center', 'FontSize', 18)

legend(b,'Medium 1','Medium 2','Medium 3', 'location','northwest');
