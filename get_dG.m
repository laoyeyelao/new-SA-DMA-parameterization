function dG_new=get_dG(temp_conv)
% The deltaG or deltaH&deltaS file
fid=fopen('deltaG_D_neutral_T0.txt','r');
fid3=fopen('deltaH_D_neutral.txt','r');
s=textscan(fid,'%f',2,'Delimiter','\n');

if ~isnumeric(s{1})
    error(['The 1st and 2nd lines in ',fn_dg,' should be the pressure and the temperature!'])
else
    pres=s{1}(1);
    temp=s{1}(2);
end
dG=0;   
% Read in the DeltaGs
% Find the number of columns
testline=fgetl(fid);
getline=fgetl(fid3);
getline=fgetl(fid3);
getline=fgetl(fid3);
while isempty(testline) || strcmp('4D',testline(1:2))~=1
    if length(dG)==4
        dG=[dG;0];
    end
    testline=fgetl(fid);
    getline=fgetl(fid3);
    posA=textscan(testline,'%s');
    dH=str2double(getline);
    dG_ori=str2double(posA{1}{2});
    dG_conv=dH-(dH-dG_ori)/temp*temp_conv;
    dG=[dG;dG_conv];
end

fclose(fid);
fclose(fid3);
dG_new=zeros(length(dG),1);
dG_new(1:5)=dG(1:5);
dG_new(6:9)=dG(6:4:18);
dG_new(10)=dG(22);
dG_new(11:14)=dG(7:4:19);
dG_new(15)=dG(23);
dG_new(16:19)=dG(8:4:20);
dG_new(20)=dG(24);
dG_new(21:24)=dG(9:4:21);
end