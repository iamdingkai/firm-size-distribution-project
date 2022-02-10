function [V2,output] = ViterateV2(w,A,pi2,Vcont2,params)

    
    [omega,beta,d,nk,nz,cf,c0,psi,k1,indexkPNoInvest2,k2,zind2,pz3]=params{:};
    %3 options, exit, continue but don't invest, continue and invest, exit,
    
    %Option 1: VoluntaryExit
    %VolExitType2 is 1 if exit by smashing and 2 if exit by hiring someone
    %to modify capital.
    temp3=zeros(nz,nk,2);
    temp3(:,:,1)=omega*k2.*(1-d)./A;
    temp3(:,:,2)=k2.*(1-d)./A-w.*c0;
    
    
    [VVolExit2,VolExit2]=max(temp3,[],3);
    VVolExit2=(1-psi).*VVolExit2+psi.*omega.*(1-d).*k2./A;

    %Option 2: No Invest
    id_linear=sub2ind([nz, nk],zind2(:),indexkPNoInvest2(:)); %[z', k']
    VcontNoInvest3=permute(repmat(reshape(Vcont2(id_linear),[nz, nk]),[1,1,nz]),[3 2 1]); %[z, k', z']
    
    VNoInvest2=pi2-cf*w+beta.*sum(VcontNoInvest3.*pz3,3);
    VNoInvest2=(1-psi).*VNoInvest2+psi.*omega.*(1-d).*k2./A;
    
%     for iz=1:nz
%         for ik=1:nk
%             VNoInvest2BF(iz,ik)=pi2(iz,ik)-cf*w+beta.*sum(Vcont2(:,indexkPNoInvest1(ik))'.*pz2(iz,:));
%         end
%     end
%     

    %Option 3: Invest or De-invest
    %This partial V is the portion of the value function that is  not
    %dependent on today's capital
    %Note that this notation is tricky, this k2 is (z,kp) NOT (z,k)
    VcontInvest3=permute(repmat(Vcont2,[1,1,nz]),[3 2 1]); %[z, k', z']
    VcontInvest2=sum(pz3.*VcontInvest3,3); %[z, k']
    
    [partofV1,indexKp1]=max(-k2./A+beta.*VcontInvest2,[],2); %[z]

    %This is back to (z,k)
    indexKpInvest2=repmat(indexKp1,1,nk); %(z,k);
    partofV2=repmat(partofV1,1,nk); %(z,k)
    VInvest2=pi2-cf*w-c0*w+k2*(1-d)./A+partofV2; %(z,k)
    VInvest2=(1-psi).*VInvest2+psi.*omega.*(1-d).*k2./A;
    %This is the amount of investment conditional on investment or
    %deinvestment happening.
    condInvestmentAmount2=(reshape(k1(indexKpInvest2),nz,nk)-k2*(1-d))./A;
    
    %Invest policy is 1 for exit, 2 for don't invest, and 3 for invest
    temp3=nan(nz,nk,3);
    temp3(:,:,1)=VVolExit2;
    temp3(:,:,2)=VNoInvest2;
    temp3(:,:,3)=VInvest2;
    [V2,ActionIndicator2]=max(temp3,[],3);

    
    %ActionIndicator2==1 means voluntary exit
    %ActionIndicator2==2 means operate but no investment
    %ActionIndicator2==3 means positive investment
    %ActionIndicator2==4 means negative investment
    ActionIndicator2((condInvestmentAmount2<0)&(ActionIndicator2==3))=4;
    %ActionIndicator2=ActionIndicator2+(investmentAmount2<0).*(ActionIndicator2==3);

    VolExit2=(ActionIndicator2==1).*VolExit2;
    %VolExit is 0 if not exiting, 1 if smashing (omega*k/A) and 2 if
    %contracting sale (k/A-c0*w)
    
    %ActionIndicator2Expanded==0 means voluntary exit and smashing
    %ActionIndicator2Expanded==1 means voluntary exit and paying fixed cost
    %ActionIndicator2Expanded==2 means operate but no investment
    %ActionIndicator2Expanded==3 means positive investment
    %ActionIndicator2Expanded==4 means negative investment
    ActionIndicator2Expanded=ActionIndicator2;
    ActionIndicator2Expanded(VolExit2==1)=0;
    %Unconditional investment amount (0 if not investing)
    InvestmentAmount2=condInvestmentAmount2.*(ActionIndicator2>=3);
    
    output={InvestmentAmount2,ActionIndicator2,ActionIndicator2Expanded,VolExit2,VVolExit2,VInvest2,VNoInvest2,indexKpInvest2,indexKp1};
end

