clc
clear
I=imread('Textile_Nom.jpg');
I=rgb2gray(I);
pixls=[128,128];
I=imresize(I,pixls);
Nom=I;
Nom=im2double(Nom);
nLevel=3;
pixs=size(Nom);
numcof=16;
wname='haar';
MWAVE=[];
for t=1:1000
    t
    Mwavecoefincontrol=[];
    U=imnoise(Nom,'gaussian',0,0.0001);
    %T=double(U)-double(Nom);
    for i=1:pixs(1)
        Detailc3=[];
        wavecoef=[];
        waveREC=[];
        x=U(i,:);
        %x=double(x);
        [C, L]=wavedec(x,nLevel,wname);
        for j=1:16
            wavecoef=[wavecoef,C(j)]; %Important coefficients
        end
        Mwavecoefincontrol=[Mwavecoefincontrol,wavecoef];% matrix of important coefficients for the whole image
    end
    MWAVE=[MWAVE;Mwavecoefincontrol];
end
MEAN=mean(MWAVE);
SQRT=std(MWAVE);
VAR=SQRT.^2;
LLL=numel(Mwavecoefincontrol);
%Phase II
RESULTS=[];
%for Z=6:10
rl=[];
changepoint=3;
RL=[];
changedev=[];
DSC=[];
Iter=100;


UCL=12.9635;

for Z=9:9
    Z
    for p=1:Iter
        p
        maxR=[];
        R=[];
        STATISTIC=0;
        counter=1;
        for H=1:changepoint
            counter=counter+1;
            Mwavecoefincontroll=[];
            MRECincontroll=[];
            U=imnoise(Nom,'gaussian',0,0.0001);
            %T=double(U)-double(Nom);
            for i=1:pixs(1)
                Detailc3=[];
                wavecoef=[];
                waveREC=[];
                x=U(i,:);
                %x=double(x);
                [C, L]=wavedec(x,nLevel,wname);
                for j=1:numcof
                    wavecoef=[wavecoef,C(j)]; %Important coefficients
                end
                Mwavecoefincontroll=[Mwavecoefincontroll,wavecoef];% matrix of important coefficients for the whole image
            end
            GLRincontroll{counter}=[];
            GLRincontroll{counter}=[Mwavecoefincontroll];
            if counter<=11
                for m=1:counter-1
                    for g=1:LLL
                        sumsum(g)=0;
                        for h=m+1:counter
                            sumsum(g)=sumsum(g)+GLRincontroll{h}(g);
                        end
                        R{counter,m}(g)=((counter-m)/(2*VAR(g)))*(((sumsum(g)/(counter-m))-MEAN(g))^2);
                        %R{counter,m}(g,l)=((counter-m)*(((sumsum(g,l)/(counter-m)))^2))/(2*VAR(g,l));
                    end
                    maxR(counter,m)=max(max(R{counter,m}));
                end
            else
                for m=counter-10:counter-1
                    for g=1:LLL
                        sumsum(g)=0;
                        for h=m+1:counter
                            sumsum(g)=sumsum(g)+GLRincontroll{h}(g);
                        end
                        R{counter,m}(g)=((counter-m)/(2*VAR(g)))*(((sumsum(g)/(counter-m))-MEAN(g))^2);
                        %R{counter,m}(g,l)=((counter-m)*(((sumsum(g,l)/(counter-m)))^2))/(2*VAR(g,l));
                    end
                    maxR(counter,m)=max(max(R{counter,m}));
                end
            end
            STATISTIC=max(maxR(counter,:));
            statistics(counter)=STATISTIC;
        end
        
        
        
        
        
        STATISTIC=0;
        while STATISTIC<UCL
            % for h=1:10
            counter=counter+1;
            Mwavecoefincontroll=[];
            MRECincontroll=[];
            IIII=Nom;
            U=imnoise(Nom,'gaussian',0,0.0001);
            for i=15:20
                for j=i:i+2
                    U(i,j)=U(i,j)+Z/255;
                end
            end
            %U=im2double(U);
            %T=double(U)-double(Nom);
            for i=1:pixs(1)
                Detailc3=[];
                wavecoef=[];
                waveREC=[];
                x=U(i,:);
                %x=double(x);
                [C, L]=wavedec(x,nLevel,wname);
                for j=1:numcof
                    wavecoef=[wavecoef,C(j)]; %Important coefficients
                end
                Mwavecoefincontroll=[Mwavecoefincontroll,wavecoef];% matrix of important coefficients for the whole image
                
            end
            GLRincontroll{counter}=[];
            GLRincontroll{counter}=[Mwavecoefincontroll];
            if counter<=11
                for m=1:counter-1
                    for g=1:LLL
                        sumsum(g)=0;
                        for h=m+1:counter
                            sumsum(g)=sumsum(g)+GLRincontroll{h}(g);
                        end
                        R{counter,m}(g)=((counter-m)/(2*VAR(g)))*(((sumsum(g)/(counter-m))-MEAN(g))^2);
                        %R{counter,m}(g,l)=((counter-m)*(((sumsum(g,l)/(counter-m)))^2))/(2*VAR(g,l));
                    end
                    maxR(counter,m)=max(R{counter,m});
                end
            else
                for m=counter-10:counter-1
                    for g=1:2048
                        sumsum(g)=0;
                        for h=m+1:counter
                            sumsum(g)=sumsum(g)+GLRincontroll{h}(g);
                        end
                        R{counter,m}(g)=((counter-m)/(2*VAR(g)))*(((sumsum(g)/(counter-m))-MEAN(g))^2);
                        %R{counter,m}(g,l)=((counter-m)*(((sumsum(g,l)/(counter-m)))^2))/(2*VAR(g,l));
                    end
                    maxR(counter,m)=max(R{counter,m});
                end
            end
            STATISTIC=max(maxR(counter,:));
            statistics(counter)=STATISTIC;
            CHANGEESTIMATE=find(maxR(counter,:)>STATISTIC-0.00000001);
            RR=[];
            % for o=1:256
            %     RR=[RR R{counter,CHANGEESTIMATE}(o,:)];
            % end
            % [a,b]=sort(RR);
            % HH=0;
        end
        RL(p)=counter-(changepoint+1);
        changedev(p)=CHANGEESTIMATE-(changepoint+1);
    end
    ARL=mean(RL);
    STDRL=std(RL);
    MEDIANCHANGE=median(changedev);
    EXPECTEDCHANGE=mean(changedev);
    STDCHANGE=std(changedev);
    RESULTS(Z,:)=[ARL STDRL MEDIANCHANGE EXPECTEDCHANGE STDCHANGE];
end