%% Figure 1
A=importdata('Modeling results\20221201_W2R_201812201901_20bin_DMA14_Mech8_Beijing.txt');
T=A.data(:,76);
p=A.data(:,77);
DMA=A.data(:,61);
CS=A.data(:,66);
SA_mono_sim=A.data(:,62)*1e6;
dH=-24.82;
dG=-13.54;
for i=1:1488
yita(i)=yita_dH_dG_mono(p(i),T(i),SA_mono_sim(i),DMA(i),CS(i),dH,dG);
end
SA=SA_mono_sim./(1-yita');
dH=-24.82;
dG=-13.54;
[J1_par,J1_ACDC,J1_dynamic]=Comparison(T,SA,DMA,CS);
sinkrate=cal_evap_rate_dH_dG(101325,T,dH,dG)+CS;
figure,
set(gcf,'Units','Centimeters','Position',[2 2 40 15]);

subplot(1,2,1)
box on
hold on 
p=scatter(J2_dynamic/1e6,J1_par/1e6,30,T);
p1=plot(10.^(-2:3),10.^(-2:3),'-k','linewidth',4);
p2=plot(10.^(-2:3),10*10.^(-2:3),'--k','linewidth',2);
plot(10.^(-2:3),10.^(-2:3)/10,'--k','linewidth',2);
p3=plot(10.^(-2:3),1.5*10.^(-2:3),'LineStyle','--','Color',[0.65 0.65 0.65],'linewidth',2);
plot(10.^(-2:3),10.^(-2:3)/2,'LineStyle','--','Color',[0.65 0.65 0.65],'linewidth',2);
legend([p,p1,p3,p2],{'\itJ\rm_{1.4}','1:1 line','\pm50% line','10:1 and 1:10 line'},'Location','southeast','Box','off');
ylabel('Parameterized \itJ\rm_{1.4} (cm^{-3}s^{-1})');
xlabel('Simulated \itJ\rm_{1.4} in KM (cm^{-3}s^{-1})');
set(gca,'yscale','log','xscale','log','xlim',[1e-2 1e3],'ylim',[1e-2 1e3],'Fontsize',20,'Fontname','Times New Roman','Ticklength',[0.025 0.025],'XMinorTick','on','TickDir','out','linewidth',1.5);
R_1=corrcoef(J2_dynamic,J1_par);
R2_1=R_1(2,1)^2;
text(0.015,400,'(a)','Fontname','Times New Roman','Fontsize',21,'FontWeight','bold');
text(0.15,400,['R^2=',num2str(R2_1)],'Fontname','Times New Roman','Fontsize',21,'FontWeight','bold');
subplot(1,2,2)
box on
hold on 
p=scatter(J1_ACDC_bk/1e6,J1_par/1e6,30,T);
p1=plot(10.^(-2:3),10.^(-2:3),'-k','linewidth',4);
p2=plot(10.^(-2:3),10*10.^(-2:3),'--k','linewidth',2);
plot(10.^(-2:3),10.^(-2:3)/10,'--k','linewidth',2);
p3=plot(10.^(-2:3),1.5*10.^(-2:3),'LineStyle','--','Color',[0.65 0.65 0.65],'linewidth',2);
plot(10.^(-2:3),10.^(-2:3)/2,'LineStyle','--','Color',[0.65 0.65 0.65],'linewidth',2);
ylabel('Parameterized \itJ\rm_{1.4} (cm^{-3}s^{-1})');
xlabel('Simulated \itJ\rm_{1.4} in CDS (cm^{-3}s^{-1})');
R_2=corrcoef(J1_ACDC_bk,J1_par);
R2_2=R_2(2,1)^2;
text(1.5e-2,400,'(b)','Fontname','Times New Roman','Fontsize',21,'FontWeight','bold');
text(0.15,400,['R^2=',num2str(R2_2)],'Fontname','Times New Roman','Fontsize',21,'FontWeight','bold');
legend([p,p1,p3,p2],{'\itJ\rm_{1.4}','1:1 line','\pm50% line','10:1 and 1:10 line'},'Location','southeast','Box','off');
c=colorbar;
c.Location='manual';
c.Position=[0.91 0.206 0.015 0.72];
set(c,'ytick',263.15:5:313.15,'yticklabel',{'-10','-5','0','5','10','15','20','25','30','35','40'},'Fontsize',16,'Fontname','Times New Roman');
ylabel(c,'Temperature (\circC)','Fontsize',21);
set(gca,'yscale','log','xscale','log','xlim',[1e-2 1e3],'ylim',[1e-2 1e3],'Fontsize',20,'Fontname','Times New Roman','Ticklength',[0.025 0.025],'XMinorTick','on','TickDir','out','linewidth',1.5);
%% Figure 2
SA=3.25e6;
T=272;
DMA=1.2;
CS=0.027;
dH=-24.819;
dG=-13.535;
fraction=10.^(-1:0.2:1)';
delta=(-20:4:20)';
J1=zeros(11,5);
J2=zeros(11,5);
dt1=zeros(11,5);
dt2=zeros(11,5);
Conc_A4B4=zeros(11,5);
Conc_A4B4_2=zeros(11,5);
varSA=SA*fraction;
J_lyy_SA=JLYYnew_dH_dG(101325,T*ones(length(varSA),1),varSA*1e6,DMA*ones(length(varSA),1),CS*ones(length(varSA),1),dH,dG);
save('.\Save_results\Para_varSA.mat','J_lyy_SA','varSA');
varDMA=DMA*fraction;
J_lyy_DMA=JLYYnew_dH_dG(101325,T*ones(length(varDMA),1),SA*ones(length(varDMA),1)*1e6,varDMA,CS*ones(length(varDMA),1),dH,dG);
yita=yita_dH_dG(101325,T*ones(length(varDMA),1),SA*ones(length(varDMA),1),varDMA,CS*ones(length(varDMA),1),dH,dG);
save('.\Save_results\Para_varDMA.mat','J_lyy_DMA','varDMA','yita');
varCS=CS*fraction;
J_lyy_CS=JLYYnew_dH_dG(101325,T*ones(length(varCS),1),SA*ones(length(varCS),1)*1e6,DMA*ones(length(varCS),1),varCS,dH,dG);
save('.\Save_results\Para_varCS.mat','J_lyy_CS','varCS');
varT=T+delta;
J_lyy_T=JLYYnew_dH_dG(101325,varT,SA*ones(length(varT),1)*1e6,DMA*ones(length(varT),1),CS*ones(length(varT),1),dH,dG);
gamma=cal_evap_rate_dH_dG(101325,varT,dH,dG);
figure,
set(gcf,'Units','Centimeters','Position',[2 2 30 22]);
subplot('Position',[0.1 0.1 0.31 0.37]);
hold on
plot(varSA,J_lyy_SA/1e6,'linewidth',3);
set(gca,'xscale','log','yscale','log','Fontsize',16,'xlim',[1e5 7e7],'Fontname','Times New Roman','Ticklength',[0.025 0.025],'XMinorTick','on','TickDir','out','linewidth',1.5);
box on
text(1.5e5,30000,'(c)','Fontname','Times New Roman','Fontsize',16,'FontWeight','bold');
xlabel('SA concentrations (cm^{-3})','Units','normalized','Position',[0.5,-0.13]);
c=ylabel('Parameterized \itJ\rm_{1.4} (cm^{-3}s^{-1})','Units','normalized','Position',[-0.12,0.5]);
subplot('Position',[0.1 0.57 0.31 0.37]);
hold on
yyaxis left
plot(varT-273.15,J_lyy_T/1e6,'linewidth',3);
ylabel('Parameterized \itJ\rm_{1.4} (cm^{-3}s^{-1})','Units','normalized','Position',[-0.12,0.5]);
set(gca,'yscale','log','Fontsize',16,'Fontname','Times New Roman','Ticklength',[0.025 0.025],'XMinorTick','on','TickDir','out','linewidth',1.5);
text(-18,500,'(a)','Fontname','Times New Roman','Fontsize',16,'FontWeight','bold');
yyaxis right
plot(varT-273.15,gamma,'linewidth',3);
ylabel('evaporation rates of \itA\rm_1\itB\rm_1 (s^{-1})','Units','normalized','Position',[1.15,0.5]);
set(gca,'yscale','log','Fontsize',16,'xlim',[-20 30],'Fontname','Times New Roman','Ticklength',[0.025 0.025],'XMinorTick','on','TickDir','out','linewidth',1.5);
 box on
xlabel('Temperature (\circC)','Units','normalized','Position',[0.5,-0.13]);

subplot('Position',[0.57 0.57 0.31 0.37]);
hold on
plot(varCS,J_lyy_CS/1e6,'linewidth',3);
set(gca,'xscale','log','yscale','log','xlim',[1e-3 3e-1],'Fontsize',16,'Fontname','Times New Roman','Ticklength',[0.025 0.025],'XMinorTick','on','TickDir','out','linewidth',1.5);
box on
xlabel('Condensation sink (s^{-1})','Units','normalized','Position',[0.5,-0.13]);
ylabel('Parameterized \itJ\rm_{1.4} (cm^{-3}s^{-1})','Units','normalized','Position',[-0.12,0.5]);
text(1.3e-3,3000,'(b)','Fontname','Times New Roman','Fontsize',16,'FontWeight','bold');
subplot('Position',[0.57 0.1 0.31 0.37]);
hold on
yyaxis left
plot(varDMA,J_lyy_DMA/1e6,'linewidth',3);
set(gca,'xscale','log','yscale','log','xlim',[1e-1 70],'Fontsize',16,'Fontname','Times New Roman','Ticklength',[0.025 0.025],'XMinorTick','on','TickDir','out','linewidth',1.5);
ylabel('Parameterized \itJ\rm_{1.4} (cm^{-3}s^{-1})','Units','normalized','Position',[-0.12,0.5]);
text(1.5e-1,500,'(d)','Fontname','Times New Roman','Fontsize',16,'FontWeight','bold');
yyaxis right
plot(varDMA,yita,'linewidth',3);
ylabel('[\itA\rm_1\itB\rm_1]/[SA_{tot}] (-)','Units','normalized','Position',[1.15,0.5]);
set(gca,'xscale','log','Fontsize',16,'Fontname','Times New Roman','Ticklength',[0.025 0.025],'XMinorTick','on','TickDir','out','linewidth',1.5,'ylim',[0 1]);
box on
xlabel('DMA concentrations (ppt)','Units','normalized','Position',[0.5,-0.13]);%%

%%
%%%#####################################################################
%%%                                                            GLOBAL                                                           
%%%#####################################################################

cases        =["DMA14_Mech8_Beijing",   "DMA17_Mech8_Beijing", "NoDMA_Mech7_Beijing", "NoDMA_Mech0_Beijing"];
cases_text   =["DMA1.4\_Mech8",   "DMA1.7\_Mech8", "NoDMA\_Mech7", "NoDMA\_Mech0"];

prefix="20221201_W2R_201812201901_20bin_";

for caca=1:length(cases)
    eval(['tmp_201812','_', char(cases(caca)),        '=', 'readtable("',         char(prefix), char(cases(caca)), '");']);
    eval(['tmp_201812','_', char(cases(caca)),'_m', '=', 'table2array(','tmp_201812','_', char(cases(caca)),  ');']);
end
clear caca prefix;

dot="'";


Readorload         = 1;
SenTest            = 0;
PSD_avg            = 1;
PSD_time           = 1;
COMPONENTS         = 0;
OTHERS             = 0;
J_analysis         = 1;
Addition           = 0;

if (SenTest)
    PSD_201812_obs=PSD_201812_obs(1:8928,:);
    DMA_201812_obs=DMA_201812_obs(1:8928,:);
    SA_201812_obs=SA_201812_obs(1:8928,:);
    J=J(1:8928,:);
    time=time(1:8928,:);
end


t_start=datenum('2018-12-01 00:00:00');t_end=datenum('2019-01-31 23:59:00'); %set start & end time!!!!!!!!!!!!!!!!!!!!!!!!!!
formatout='yyyy-mm-dd HH:MM:SS';
interval=datenum(hours(120));
S=timerange(datestr(t_start,formatout), datestr(t_end,formatout));
Diameter_sim_14=Diameter_sim;
Diameter_sim_14(1)=1.4;Diameter_sim_14(2)=1.7;
Diameter_sim=Diameter_sim_14;

if (Readorload) %chang here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


time_PSD_obs=t_start:datenum(minutes(5)):t_end; 
time_CompOt_obs=t_start:datenum(hours(1)):t_end;
time_DMASA_obs=t_start:datenum(minutes(5)):t_end;
time_sim =t_start:datenum(hours(1)):t_end;

end
%%%#####################################################################
%%%                                                            PSD_avg                                                           
%%%#####################################################################
if(PSD_avg)
for caca=1:length(cases)
     eval(['PSD_sim_', char(cases(caca)), '_AN', '=tmp_201812_',   char(cases(caca)), '_m(:,  1:20);']);
     eval(['PSD_sim_', char(cases(caca)), '_AC', '=tmp_201812_',   char(cases(caca)), '_m(:,21:40);']);
     eval(['PSD_sim_', char(cases(caca)), '_0',  '=tmp_201812_',   char(cases(caca)), '_m(:,41:60);']);
	 eval(['PSD_sim_', char(cases(caca)), '_F',  '=', '(PSD_sim_', char(cases(caca)), '_AN  +', 'PSD_sim_', char(cases(caca)), '_AC +', 'PSD_sim_', char(cases(caca)), '_0)./3;',]);
end
clear caca;

figure(1);
PSD_201812_obs_com=PSD_201812_obs(1:round(length(PSD_201812_obs)/length(time_sim)):length(PSD_201812_obs),:);
cri=round(length(Diameter_obs)*0.4);% 0.4 is uesd now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
row_201812=find(sum(isnan(PSD_201812_obs_com),2)<cri);
PSD_201812_obs_avg=mean(PSD_201812_obs_com(row_201812,:),1,"omitnan");

fPSD_avg=tiledlayout(1,1,'TileSpacing','compact');
fPSD_avg.Position=[0.15 0.15 0.55 0.8];

nexttile
loglog(Diameter_obs, PSD_201812_obs_avg,'LineStyle','-', 'LineWidth',8);
hold on;
for rr=1:length(cases)
     case_now=['PSD_sim_', char(cases(rr)), '_F'];
     disp(['PSD_avg: '  case_now]);
     eval([case_now, '_avg=', 'mean(', char(case_now), '(row_201812,:),1,"omitnan");' ]);
     if (rr == 1)
         loglog(Diameter_sim, eval([case_now, '_avg']),'LineStyle','-', 'Color', [0.85 0.33 0.10], 'LineWidth',8);
     elseif (rr==2)
         loglog(Diameter_sim, eval([case_now, '_avg']),'LineStyle',':', 'Color', [0.85 0.33 0.10],'LineWidth',8);
     elseif(rr==3)
         loglog(Diameter_sim, eval([case_now, '_avg']),'LineStyle','-',  'Color', [1 0.65 0], 'LineWidth',8);
      elseif(rr==4)
         loglog(Diameter_sim, eval([case_now, '_avg']),'LineStyle','-',  'Color', [0.55 0.54 0.54], 'LineWidth',8);        
     end 
     ylim([1 100000]);xlim([1 3000])
     hold on;
end

     ylim([1 100000]);
     set(gca,'tickdir', 'out',  'TickLength',[0.02 0.025],  'YTick', [1 10 100 1000 10000 100000]);
     set(gca,'FontSize',22,'Fontname','Times New Roman');

box on;set(gca,'LineWidth',2);
legend(["observation"',cases_text],'box','off')
clear rr case_now_1 case_now;

xlabel('{\itD}_p (nm)','FontSize',26,'Fontname','Times New Roman');
ylabel( 'd{\itN}/dlog{\itD}_p (# cm^{-3})',FontSize=26,Fontname='Times New Roman');
set(gcf,'outerposition',get(0,'screensize'));
% exportgraphics(fPSD_avg, 'PNSD_201812_avg.png','resolution', 900);
end

%%%#####################################################################
%%%                                                            PSD_time                                                           
%%%#####################################################################
if(PSD_time)
figure(2);
fPSD201812=tiledlayout(length(cases)+1,1,TileSpacing="compact");
fPSD201812.Position=[0.1 0.1 0.55 0.85];
for rr=1:length(cases)+1
     nexttile
     if(rr==1)
         pcolor(time_PSD_obs,Diameter_obs,log10(PSD_201812_obs'));colormap(jet);caxis([0 6.5]);grid off;shading flat;
         set(gca,'XTick',t_start:interval:t_end);datetick('x',6,'keeplimits','keepticks');
         set(gca,'xticklabels',[]);
         set(gca,'YScale','log');
         ylim([1 10000]);
         set(gca,'FontSize',16,'Fontname','Times New Roman');
         title('\bf observation','FontSize',16,'Fontname','Times New Roman');
         set(gca,'TickDir','out','XTick',[]);
         box off;
     else
         case_now=['PSD_sim_', char(cases(rr-1)), '_F'];
         case_now_1=char(cases_text(rr-1));
         disp(['PSD_time: '  case_now]);
         pcolor(time_sim,Diameter_sim,log10(eval([char(case_now),  char(dot)])));colormap(jet);caxis([0 6.5]);grid off;shading flat;
         set(gca,'xticklabels',[]);
         set(gca,'YScale','log');
         ylim([1 10000]);
         set(gca,'FontSize',16,'Fontname','Times New Roman');
         title(['\bf ' case_now_1],'FontSize',16,'Fontname','Times New Roman');
         set(gca,'TickDir','out','XTick',[]);
         box off;
     end
     if(rr==length(cases)+1)
         set(gca,'XTick',t_start:interval:t_end,'FontSize', 16,'Fontname','Times New Roman');datetick('x',6,'keeplimits','keepticks');
         title(['\bf' case_now_1],'FontSize',16,'Fontname','Times New Roman');
     end
     set(gca,'linewidth',1);
end
clear rr case_now case_now_1;



xlabel(fPSD201812,'Date','FontSize',20,'Fontname','Times New Roman');
ylabel(fPSD201812,'{\itD}_p (nm)','FontSize',20,'Fontname','Times New Roman');

cb1=colorbar;
ylabel(cb1,'d{\itN}/dlog{\itD}_p (# cm^{-3})',FontSize=16,Fontname='Times New Roman');
set(cb1,'Ticks', [2, 4, 6], 'TickLabels', {'10^2', '10^4', '10^6'});
cb1.Layout.Tile='east';
 



set(gcf,'outerposition',get(0,'screensize'));
% exportgraphics(fPSD201812,'PNSD_201812.png','resolution', 900);
end


%%%#####################################################################
%%%                                                            OTHERS                                                         
%%%#####################################################################
if(OTHERS)
ot=[ "DMA", "SA"];
for caca=1:length(cases)
     eval(['DMA_sim_201812',    char(cases(caca)), '=tmp_201812_', char(cases(caca)), '_m(:,61);']);
     eval(['SA_sim_201812',       char(cases(caca)), '=tmp_201812_', char(cases(caca)), '_m(:,62);']);
end
clear caca;

figure(3);

fot201812=tiledlayout(2,1,TileSpacing="compact");
fot201812.Position=[0.12 0.12 0.8 0.8];
for zz=1:length(ot)
    nexttile;
    case_now=[char(ot(zz)) '_201812_obs'];
    eval(['scatter(time_DMASA_obs, ', char(case_now), ', 8, LineWidth=1)']);hold on;

for zzz=1:1%length(cases)
    name=[ char(ot(zz)),'_sim_201812', char(cases(zzz))];
    eval(['plot(time_sim, ', char(name), ',LineWidth=3)']);hold on;
    xlim([t_start t_end]);

if (zz==4) % here may be chaged!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   legend_cases=["observation" cases_text];set(gca,'XTick',t_start:interval:t_end,'xticklabels',[]);set(gca,'FontSize',30,'Fontname','Times New Roman');
   yl=ylabel('\bf PM_{25} (ug cm^{-3})',FontSize=14,Fontname='Times New Roman');legend(legend_cases(:), 'Orientation','horizontal','Location','northoutside');
elseif(zz==1)
   set(gca,'XTick',t_start:interval:t_end,'xticklabels',[]);set(gca,'FontSize',22,'Fontname','Times New Roman');yl=ylabel('DMA (ppt)',FontSize=26);
%    text('string', '\bf (a)','Units', 'normalized', 'position',[0.05 0], 'FontSize',26,'Fontname','Times New Roman');
elseif(zz==2)
   set(gca,'XTick',t_start:interval:t_end,'xticklabels',[]);set(gca,'FontSize',22,'Fontname','Times New Roman');yl=ylabel('SA (# cm^{-3})',FontSize=26);datetick('x',6,'keeplimits','keepticks');
%    text('string', '\bf (b)','Units', 'normalized', 'position',[0.05 0], 'FontSize',26,'Fontname','Times New Roman');
else
   set(gca,'XTick',t_start:interval:t_end);datetick('x',6,'keeplimits','keepticks');set(gca,'FontSize',16,'Fontname','Times New Roman');yl=ylabel('\bf SO_2 (ug cm^{-3})',FontSize=14,Fontname='Times New Roman');
end
 
clear yl mean_1 mean_2 mean_obs mean_sim mean_sim_1;
set(gca,'tickdir', 'out',  'TickLength',[0.01 0.025]);box on;set(gca,'LineWidth',2);

if(zz==1)
    legend('observation','DMA1.4\_Mech8','box','off');
end

end
end



clear zz zzz;
xlabel('Date','FontSize',26,'Fontname','Times New Roman');
set(gcf,'outerposition',get(0,'screensize'));
set(gca,'LooseInset',get(gca,'TightInset'));
end

%%%#####################################################################
%%%                                                            J ANALYSIS                                                          
%%%#####################################################################
% if(J_analysis)

J_sim_obs = tmp_201812_DMA14_Mech8_Beijing_m(:,75);
J_1h=J(1:12:end);
% for xxx=1:length(J_1h)
%     if (isnan(J_1h(xxx)))
%         J_sim_obs(xxx) = nan;
%     end
% end

for xxx=1:size(PSD_201812_obs_com,1)
    if (sum(isnan(PSD_201812_obs_com(xxx,:))) >= 0.4*size(PSD_201812_obs_com,2) )
        J_sim_obs(xxx)=nan;
    end
end


if(J_analysis)
Mech8=tmp_201812_DMA14_Mech8_Beijing_m(:,[68:74,67]);
Mech7=tmp_201812_DMA14_Mech8_Beijing_m(:,68:74);
Mech7(:,8)=nan;

figure(4);

fj=tiledlayout(2,2,TileSpacing="compact");
%  fj=tiledlayout('flow');
% fj.Position=[0.10 0.12 0.7 0.8];
nexttile(1, [1 2])
% a1=area(time_sim,Mech8);a1(1).FaceColor=[0 0 0];a1(1).FaceColor=[0 0 0];a1(2).FaceColor=[0.55 0.53 0.44];a1(8).FaceColor=[0.85 0.33 0.10];hold on;
scatter(time, (J./1e6),22,[0.0 0.45 0.74],'LineWidth',1);hold on;
% plot(time_sim,tmp_201812_DMA14_Mech8_Beijing_m(:,75),LineWidth=2,Color= [0.85 0.33 0.10]);
plot(time_sim,J_sim_obs,LineWidth=2,Color= [0.85 0.33 0.10]);

xlim([t_start t_end]); set(gca,'FontSize',22,'Fontname','Times New Roman');
ylabel('{\it{J}} (#/cm^{-3} s^{-1})', 'FontSize',26,'Fontname','Times New Roman');
set(gca,'tickdir', 'out',  'TickLength',[0.01 0.025]);box on;set(gca,'LineWidth',2);
 set(gca,'XTick',t_start:interval:t_end);datetick('x',6,'keeplimits','keepticks');xlabel('Date','FontSize',26,'Fontname','Times New Roman');
legend('observation','DMA1.4\_Mech8','box','off');

% text('string', '\bf (a)','Units', 'normalized', 'position',[0.05 0], 'FontSize',26,'Fontname','Times New Roman');


nexttile(3);

C = [068 004 090
       065 062 133
       048 104 141
       031 146 139
       053 183 119
       145 213 066
       248 230 032]/255;

C1= [0.5 0 0
        [96 96 96]/255];

pp1=pie3(sum([Mech8(:,8) sum(Mech8(:,1:7),2)]));colormap(C1);


legend('DMA+sulfuric acid','others','box','off',Location='eastoutside');
set(gca,'FontSize',22,'Fontname','Times New Roman');
% le.Layout.Tile='east';
nt2=nexttile(4);
pp2=pie3(sum(Mech8(:,1:7),1));colormap(nt2,C);
legend('binary neutral', 'binary ion-induced', 'ternary neutral','ternary ion-indued','organic neutral','organic ion-induced','organic+sulfuric acid','box','off',Location='eastoutside');
set(gca,'FontSize',22,'Fontname','Times New Roman');


end











