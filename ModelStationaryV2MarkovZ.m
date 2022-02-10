
clc;
clear variables;


c0vec=3; %4
cEvec=10; %4
cfvec=1; %3
psivec=0.05; %1
phivec=0.03; %2
Avec=0.1:0.1:1;

[c0mat, cEmat, cfmat, psimat, phimat, Amat]=ndgrid(c0vec,cEvec,cfvec,psivec,phivec,Avec);
c0long=c0mat(:);
cElong=cEmat(:);
cflong=cfmat(:);
psilong=psimat(:);
philong=phimat(:);
Along=Amat(:);
iterations=length(c0long);
iterMax=100;

for loopCount=1:iterations
    loopCount
    tic
    %try
        flag1=0;
        flag2=0;
        c0=c0long(loopCount); %Fixed labor cost of investment
        cE=cElong(loopCount); %Entry fixed cost
        cf=cflong(loopCount); %Continuation fixed cost
        psi=psilong(loopCount); %psi is exogenous probability of firm exit.
        phi=philong(loopCount); %Scalar on disultility of labor
        A=Along(loopCount); %Scalar the rate that final good can be turned into capital good.
        
        display(strcat('c0=',num2str(c0),' cE=',num2str(cE),' cf=',num2str(cf),' psi=',num2str(psi),' phi=',num2str(phi)));

        stepPercent=0.8; %How fast to close in on w vector
        %Consumer utility variables.
        nu=0.5; %Frisch elasticity of labor paramater
        sigma=2; %intertermporal elasticity of substitution
        uc=@(C)C.^(-sigma).*(C>0)+(C<=0).*1000000000.*(-C);
        un=@(N)phi*N.^(nu);
        eta=0.95; %persistence of efficiency of investment parameter
        omega=0.0; %Fraction of specific capital that can be salvaged
        minMeasure=0.000; %Minimum measure of firms to graph

        %alpha is the decreasing returns to scale parameter of the production
        %function.
        alpha=0.4;
        xi=0.85; %Decreasing return to scale paramter
        %Discount parameter
        beta=0.96;

        %depreciation
        d=0.05;

        %Number of grid points in k
        nk=800;
        kL=0.001; %Smallest Capital
        kEMin=kL;
        % kE=1.5; %Entry capital fraction of kU

        %Markov pattern for z.
        nz=50;
        
%         REMOVE THIS WHEN WE'RE DONE
%         pz2=1;
%         z1=1;
%         logz1=log(z1);
%         pzstationary1=1;
%         pzstationary2=pzstationary1.*ones(1,nk);
%         REMOVE THIS WHEN WE'RE DONE
        
% THIS IS THE TAUCHEN CODE PUT IT BACK IN WHEN WE'RE DONE
        mulogz=0;
        rhologz=.8; %persistance for log z
        sigmalogz=0.3;
        mlogz=3.5;


        [logz1,pz2] = tauchen(nz,mulogz,rhologz,sigmalogz,mlogz);
        z1=exp(logz1);
        pz3=permute(repmat(pz2,[1,1,nk]),[1 3 2]); %[z, k', z']
        pzstationary1=(ones(1,nz)./nz*pz2^200)';
        [pzstationary2,~]=ndgrid(pzstationary1,1:nk);
% THIS IS THE TAUCHEN CODE PUT IT BACK IN WHEN WE'RE DONE

        %If persistence is 1 use this code
        % [logz1,pz2] = tauchen(nz,mulogz,0,sigmalogz,mlogz);
        % z1=exp(logz1);
        % pzstationary1=pz2(1,:)';
        % pz2=eye(nz);
        % pz3=permute(repmat(pz2,[1,1,nk]),[1 3 2]); %[z, k', z']
        % [pzstationary2,~]=ndgrid(pzstationary1,1:nk);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1212)
        plot(logz1,pzstationary1)

        % %Number of grid points in z
        % nz=100;
        % zL=1;
        % z1=zL*ones(1,nz);
        % gamma=3;
        % %Pareto PDF: gamma*zmin^gamma/(z^(1+gamma)
        % %Pareta CDF: 1-(zL/z)^gamma;
        % %i=1-(zL/z)^gamma (Solve for z)
        % %z=zL/(1-i)^(1/gamma)
        % for i=1:(nz-1)
        %     z1(i+1)=zL/(1-i/nz)^(1/gamma);
        % end
        % fz1=(1/nz).*ones(1,nz);





        %New  Version of investment to capital good



        % nt=max(1000,-log(c1zero/c1bar)/log(eta))

        %Strategy: Guess a T.w, solve value function, iterate and verify T.w is
        %supported.
        % wMax=20;
        % wMin=0.1;
        % w=(wMax+wMin)/2;
        w=0.959858716519416;
        wMax=40;
        wMin=0.2;
        VE=2*cE*w;
        indexkE=nan;
        %[fz2,~]=ndgrid(fz1,1:nk);
        kU=((z1(nz)*xi)^(1/(1-xi))*(alpha*A/d)^((1-xi*(1-alpha))/(1-xi))*((1-alpha)/w)^(((1-alpha)*xi)/(1-xi)))/5;
        iterCount=0;
        while max(abs(VE-cE*w))>0.000001 && iterCount<iterMax
            iterCount=iterCount+1;
            kUUpdate=((z1(nz)*xi)^(1/(1-xi))*(alpha*A/d)^((1-xi*(1-alpha))/(1-xi))*((1-alpha)/w)^(((1-alpha)*xi)/(1-xi)))/5;
            kU=(abs(kU-kUUpdate)>10000)*(kUUpdate-kU)+kU;
            logk1=linspace(log(kL),log(kU),nk)';
            k1=exp(logk1);

            display(strcat('     entryCond=',num2str(VE-cE*w),' wMin=',num2str(wMin),' wMax=',num2str(wMax),' kU=',num2str(kU)));
            
        %     logk1part1=linspace(log(kL),log(kU),nk/2)';
        %     k1part1=exp(logk1part1);
        %     k1part2=linspace(kL,kU,nk/2)';
        %     k1=sort(unique([k1part1; k1part2]));
        %     nktemp=length(k1);
        %     k1(nktemp+1:nk)=linspace(kU+1e-5,kU+1,nk-nktemp);

            %Grid of k
            kPNoInvest1=k1.*(1-d);
            indexkPNoInvest1=knnsearch(k1,kPNoInvest1);
            indexkPNoInvest2=repmat(indexkPNoInvest1',nz,1);
            %Calculate stationary equilibrium
            %Make matrices for stationary Equilibrium
            [z2, k2]=ndgrid(z1,k1);
            [zind2, kind2]=ndgrid(1:nz,1:nk);
            Vcont2=ones(nz,nk);
            V2=zeros(nz,nk);

            %nPolicy is just the labor used towards production, not labor hired for
            %fixed costs
            nPolicy2=((z2*xi*(1-alpha).*k2.^(alpha*xi))./w).^(1/(1-xi+alpha*xi));
            q2=z2.*(k2.^alpha.*nPolicy2.^(1-alpha)).^xi;
            pi2=q2-w.*nPolicy2;
            params={omega,beta,d,nk,nz,cf,c0,psi,k1,indexkPNoInvest2,k2,zind2,pz3};

            %% solve for stationary equilibrium given w
            while(max(max(abs(Vcont2-V2)))>0.0001)
                Vcont2=V2;
                [V2,output]=ViterateV2(w,A,pi2,Vcont2,params);
            end
            [investmentAmount2,ActionIndicator2,ActionIndicator2Expanded,VolExit2,VVolExit2,VInvest2,VNoInvest2,indexKpInvest2,indexKp1]=output{:};
            indexKp2=(ActionIndicator2==2).*indexkPNoInvest2+(ActionIndicator2>=3).*indexKpInvest2;
            indexKp2(VolExit2>=1)=nan; %In reality it's zero, but we don't want to mess this up with an actual 0 so it will error if we try to use it.

            %VE is entry value, note: involuntary exit may happen at the beginning
            %of next period but it's already included in V2
            [VE,indexkE]=max(-cE*w-c0*w-k1./A+beta*sum(V2.*pzstationary2,1)'); 

            %This code makes a mandotory minimum starting capital level.
        %     if k1(indexkE)<kEMin
        %         indexkE=find(k1>kEMin,1,'first');
        %         temp=beta*sum(V2.*fz2,1)';
        %         VE=-cE*w-c0*w-k1(indexkE)./A+temp(indexkE);
        %     end

            wMin=(VE>cE*w)*(stepPercent*w+(1-stepPercent)*wMin)+(VE<=cE*w)*wMin;
            wMax=(VE<cE*w)*(stepPercent*w+(1-stepPercent)*wMax)+(VE>=cE*w)*wMax;
            w=(wMin+wMax)/2;
        end

        if iterCount==iterMax
            flag1=1;
        end


        % figure(3);clf;
        % hold on
        % plot(k1,VE2(:,13),'bo-',k1,VE2(:,14),'rx-');
        % xlim([0 k1(300)]);
        % plot(k1(300),VE2(300,13),'k+')
        % plot(k1(241),VE2(241,14),'k+')
        % hold off

        figure(1);clf;
        subplot(2,3,1);
        subplot(2,3,1);
        surf(log(z2(:,:)),log(k2(:,:)),VInvest2(:,:))
        xlabel('log(z)')
        ylabel('log(k)')
        zlabel('VInvest2')
        subplot(2,3,2);
        surf(log(z2(:,:)),log(k2(:,:)),VNoInvest2(:,:))
        xlabel('log(z)')
        ylabel('log(k)')
        zlabel('VNoInvest2')
        subplot(2,3,3);
        surf(log(z2(:,:)),k2(:,:),V2(:,:))
        xlabel('log(z)')
        ylabel('log(k)')
        zlabel('V2')
        subplot(2,3,4);
        surf(log(z2(:,:)),log(k2(:,:)),ActionIndicator2Expanded(:,:))
        xlabel('log(z)')
        ylabel('log(k)')
        zlabel('InvestPolicy')
        subplot(2,3,5);
        surf(log(z2(:,:)),log(k2(:,:)),log(k1(indexKpInvest2(:,:))))
        xlabel('log(z)')
        ylabel('log(k)')
        zlabel('log(kp)|invest')
        subplot(2,3,6);
        surf(log(z2(:,:)),log(k2(:,:)),nPolicy2(:,:))
        xlabel('log(z)')
        ylabel('log(k)')
        zlabel('labor')

        zIndex=round(nz/2);
        figure(2);clf;
        subplot(2,2,1);
        hold on
        plot(log(k1),log(VInvest2(zIndex,:)))
        plot(log(k1),log(VNoInvest2(zIndex,:)))
        plot(log(k1),log(VVolExit2(zIndex,:)))
        plot(log(k1),log(V2(zIndex,:)))
        legend('Invest','No Invest','Exit','Value')
        xlabel('k')
        ylabel('V')
        subplot(2,2,2);
        plot(log(k1),ActionIndicator2Expanded(zIndex,:))
        xlabel('k')
        ylabel('InvestPolicy')
        subplot(2,2,3);
        plot(log(k1),log(k1(indexKpInvest2(zIndex,:))))
        xlabel('k')
        ylabel('log(kp)|invest')
        subplot(2,2,4);
        plot(log(k1),investmentAmount2(zIndex,:))
        xlabel('k')
        ylabel('investAmount')

        figure(10);clf;
        surf(log(z2(:,1:nk/40:nk)),log(k2(:,1:nk/40:nk)),ActionIndicator2Expanded(:,1:nk/40:nk))
        xlabel('log(z)')
        ylabel('log(k)')
        zlabel('InvestPolicy')
        colorbar
        title('0 - exit do not sell capital; 1 - exit pay fixed cost and sell capital; 2 - operate no invest; 3 - operate invest; 4 - operate divest')
        hold on
        plot3(log(z1),log(k1(indexKpInvest2(:,1))),4*ones(nz,1),'bo-');
        xlabel('log(z)')
        ylabel('log(k)')
        zlabel('log(kp)|invest')

        figure(11);clf;
        surf(log(z2(:,:)),log(k2(:,:)),log(k1(indexKpInvest2(:,:))))
        xlabel('log(z)')
        ylabel('log(k)')
        zlabel('log(kp)|invest')
        %     surf(log(z2),log(k2),log(k1(indexKp2)))
        %     xlabel('log(z)')
        %     ylabel('log(k)')
        %     zlabel('log(kp)')

        %Calculate measureof firms (need a sequence of M0(t)).
        %Approach, assume we start in a steady state of c1 to get a steady state
        %M0(1). And then iterate from the.

        %% Calculate the measure of firms in the statioanry distribution
        continueLinearIndex=find(~isnan(indexKp2));
        [zindex,kindex]=ind2sub(size(indexKp2),continueLinearIndex);

        %Make a list of kp used
        temp=unique(indexKp2(:));
        kpList=temp(~isnan(temp));
        izList2=cell(nz,length(kpList)); %zp, kpList
        temp2=cell(nz,length(kpList));
        for izp=1:nz
            for j=1:length(kpList)
                ikp=kpList(j);
                temp2{izp,j}=find(indexKp2==ikp); % give all combinations of (iz,ik) that goes to ikp
                [izList2{izp,j},~]=ind2sub([nz,nk],temp2{izp,j}); 
            end
        end

        input={indexkPNoInvest1;kpList;izList2;temp2;VolExit2;pzstationary1;continueLinearIndex;zindex;kindex;indexKp2;ActionIndicator2;investmentAmount2;nPolicy2;q2;k1;k2;w;uc;un;nz;nk;psi;cf;c0;cE;omega;d;A;indexkE;pz2;indexKp1};
        [M0,fval]=fsolve(@(M0)InitialDistributionM0SolverV2(M0,input),1,optimoptions('fsolve','MaxIterations',100,'display','iter','tolfun',1e-7));
        [~,I]=InitialDistributionM0SolverV2(M0,input);

        if abs(fval)>=1e-7
            flag2=abs(fval);
        end

        logz=log(z2(:));logk=log(k2(:));
        MColumn=I.M2(:);
        figure(2);clf
        hold on
        scatter3(logz(MColumn>minMeasure),logk(MColumn>minMeasure),MColumn(MColumn>minMeasure))
        % patch(log([z1(1) z1(nz) z1(nz) z1(1)]), log([k1(indexkE) k1(indexkE) k1(indexkE) k1(indexkE)]), [0 0 max(max(M2)) max(max(M2))], 'r', 'FaceAlpha', 0.5)
        xlabel('log(z)')
        ylabel('log(k)')
        zlabel('M')
        view([65 35])

        %Remember: ActionPolicy2==4 means divesting.
        ActionIndicatorColumn=ActionIndicator2Expanded(:);
        figure(6);clf
        scatter3(logz(MColumn>minMeasure),logk(MColumn>minMeasure),ActionIndicatorColumn(MColumn>minMeasure))
        % patch(log([z1(1) z1(nz) z1(nz) z1(1)]), log([k1(indexkE) k1(indexkE) k1(indexkE) k1(indexkE)]), [0 0 max(max(M2)) max(max(M2))], 'r', 'FaceAlpha', 0.5)
        xlabel('log(z)')
        ylabel('log(k)')
        zlabel('ActionPolicy')
        view([65 35])

        M2=I.M2;
        %Calculate moments
        %Total Exit Rate
        %psi is exogenous exit rate.

        %c0 cost we're thinking of as consultants, they're not part of the firm for
        %firm size distribution moment matching purposes.
        fLaborSize2=nPolicy2+cf;
        laborGroupCutoffs=[0 4 9 19 99 499 999 2499 4999 9999 inf];
        fExitRatebySize1=nan(length(laborGroupCutoffs)-1,1);
        fMeasurebySize1=nan(length(laborGroupCutoffs)-1,1);
        for i=1:length(laborGroupCutoffs)-1
            fMeasurebySize1(i)=sum(sum((fLaborSize2>=laborGroupCutoffs(i)).*(fLaborSize2<laborGroupCutoffs(i+1)).*M2));
            fExitRatebySize1(i)=psi+(1-psi)*sum(sum((fLaborSize2>=laborGroupCutoffs(i)).*(fLaborSize2<laborGroupCutoffs(i+1)).*M2.*(VolExit2>=1)))/fMeasurebySize1(i);
        end
        fExitRatebySize1
        measureOfFirms=sum(sum(M2))
        operatingMeasureofFirms=sum(sum((1-psi).*M2.*(VolExit2==0)))
        fExitRate=psi+(1-psi)*sum(sum((VolExit2>=1).*M2))/measureOfFirms

        averageFirmSize=sum(sum(fLaborSize2.*M2))./measureOfFirms
        positiveSpikeShare=(sum(sum((ActionIndicator2==3).*M2))+M0)/measureOfFirms
        negativeSpikeShare=sum(sum((ActionIndicator2==4).*M2))/measureOfFirms
        momAndPopShare=sum(sum((fLaborSize2<laborGroupCutoffs(2)).*(ActionIndicator2~=3).*M2))/measureOfFirms
        laborVariance=sum(sum((fLaborSize2-averageFirmSize).^2.*M2))/measureOfFirms
        totalEmployment=sum(sum(fLaborSize2.*M2));
        EmploymentShares1=nan(length(laborGroupCutoffs)-1,1);
        EmploymentShares1(1)=sum(sum((fLaborSize2<laborGroupCutoffs(2)).*M2.*fLaborSize2))./totalEmployment;
        EmploymentShares1(10)=sum(sum((fLaborSize2>laborGroupCutoffs(10)).*M2.*fLaborSize2))./totalEmployment;

        for i=2:length(laborGroupCutoffs)-2
            EmploymentShares1(i)=sum(sum((fLaborSize2>=laborGroupCutoffs(i)).*(fLaborSize2<laborGroupCutoffs(i+1)).*M2.*fLaborSize2))./totalEmployment;
        end

        fExitRateFit=abs((fExitRate/0.0818)-1);
        averageFirmSizeFit=abs((averageFirmSize/24.5)-1);
        positiveSpikeShareFit=abs((positiveSpikeShare/0.186)-1);
        momAndPopShareFit=abs((momAndPopShare/0.584)-1);
        momAndPopEmploymentShareFit=abs((EmploymentShares1(1)/0.052258)-1);
        fittingDistance=(sqrt(fExitRateFit)+sqrt(averageFirmSizeFit)+sqrt(positiveSpikeShareFit)+sqrt(momAndPopShareFit)+sqrt(momAndPopEmploymentShareFit))^2;

%         params2=[omega beta d nk nz cE cf c0 psi xi nu sigma phi eta alpha A mulogz rhologz sigmalogz mlogz ...
%             fExitRatebySize1' measureOfFirms operatingMeasureofFirms fExitRate averageFirmSize positiveSpikeShare negativeSpikeShare momAndPopShare ...
%             laborVariance totalEmployment EmploymentShares1' w fittingDistance flag1 flag2];
%         xlsappend('Calibrationresults.xlsx',params2,'Results')

        empSharesData=[0.05227559	0.052258215	0.065548037	0.163709299	0.140307943	0.054294827	0.071735237	0.055675666	0.056545401	0.287649786];
        empSharesCumulativeData=[0.05227559	0.104533805	0.170081842	0.333791141	0.474099084	0.528393911	0.600129148	0.655804813	0.712350214	1];

        save(strcat('TwoZAis',num2str(floor(A*10)),'.mat'));
        figure(30);
        clf;
        hold on;
        X = categorical({'0-4','5-9','10-19','20-99','100-499','500-999','1000-2499','2500-4999','5000-9999','10000+'});
        X = reordercats(X,{'0-4','5-9','10-19','20-99','100-499','500-999','1000-2499','2500-4999','5000-9999','10000+'});
        Y=[cumsum(empSharesData') cumsum(EmploymentShares1)];
        bar(X,Y)

        figure(31);
        clf;
        hold on;
        X=categorical({'0-100','100+'})
        temp1=[sum(empSharesData(1:4)) sum(empSharesData(5:end))];
        temp2=[sum(EmploymentShares1(1:4)) sum(EmploymentShares1(5:end))];
        Y=[temp1' temp2'];
        bar(X,Y)
        hold off;

%     catch
%         c0=c0long(loopCount); %Fixed labor cost of investment
%         cE=cElong(loopCount); %Entry fixed cost
%         cf=cflong(loopCount); %Continuation fixed cost
%         psi=psilong(loopCount); %psi is exogenous probability of firm exit.
%         phi=philong(loopCount); %Scalar on disultility of labor
% 
%         stepPercent=0.8; %How fast to close in on w vector
%         %Consumer utility variables.
%         nu=0.5; %Frisch elasticity of labor paramater
%         sigma=2; %intertermporal elasticity of substitution
%         uc=@(C)C.^(-sigma).*(C>0)+(C<=0).*1000000000.*(-C);
%         un=@(N)phi*N.^(nu);
%         eta=0.95; %persistence of efficiency of investment parameter
%         omega=0.0; %Fraction of specific capital that can be salvaged
%         minMeasure=0.000; %Minimum measure of firms to graph
% 
%         %alpha is the decreasing returns to scale parameter of the production
%         %function.
%         alpha=0.4;
%         xi=0.85; %Decreasing return to scale paramter
%         %Discount parameter
%         beta=0.96;
% 
%         %depreciation
%         d=0.05;
% 
%         %Number of grid points in k
%         nk=800;
%         kL=0.001; %Smallest Capital
%         kEMin=kL;
%         % kE=1.5; %Entry capital fraction of kU
% 
%         %Markov pattern for z.
%         nz=50;
%         mulogz=0;
%         rhologz=.8; %persistance for log z
%         sigmalogz=0.3;
%         mlogz=3.5;
%         A=1; %The rate that final good can be turned into capital good.
%         params2=[omega beta d nk nz cE cf c0 psi xi nu sigma phi eta alpha A mulogz rhologz sigmalogz mlogz];
%         xlsappend('Calibrationresults.xlsx',params2,'Results')
%     end
    toc
end
