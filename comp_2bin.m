clear all
close all
clc

N=1;
sz=1;

s1=0.3; 
s2=0.4;
% b1=0.5; 
% b2=b1;

% mumin=max([-b1/s1,-b2/s2]); 
% mumax=5; 
% mu=linspace(mumin,mumax,sz);

% mu=1;
n1=1;
m1=3;
n2=2;
m2=3;

In1=NaN(N,sz); Im1=NaN(N,sz); In2=NaN(N,sz); Im2=NaN(N,sz);
% B1_f=NaN(N,sz); B2_f=NaN(N,sz); MU_f=NaN(N,sz); Lmin_f=NaN(N,sz);
MU_fm=NaN(N,sz); B1_fm=NaN(N,sz); B2_fm=NaN(N,sz); Lmin_fm=NaN(N,sz); comp_ffm=NaN(N,sz);
% check_ffm=NaN(1,sz);

%%
for y=1:sz
%     In1(:,y)=poissrnd(b1,N,1);
%     Im1(:,y)=poissrnd(b1+mu(y)*s1,N,1);
%     In2(:,y)=poissrnd(b2,N,1);
%     Im2(:,y)=poissrnd(b2+mu(y)*s2,N,1);
    In1=n1*ones(N,sz); Im1=m1*ones(N,sz); In2=n2*ones(N,sz); Im2=m2*ones(N,sz);

    %% Fit b1 b2 and mu
%     x0=[1,1,10]; lb=[-Inf,0,0]; options=optimoptions('fmincon','Display','off');
%     for i=1:N
%         n1=In1(i,y); m1=Im1(i,y); n2=In2(i,y); m2=Im2(i,y); %s1=Is1(i); s2=Is2(i);
%         A=[-s1,-1,0;-s2,0,-1;0,0,0];
%         [X,fval]=fmincon(@(C)-log(C(2)^n1/factorial(n1)*exp(-C(2))*(C(2)+C(1)*s1)^m1/factorial(m1)*exp(-C(2)-C(1)*s1)) ...
%                              -log(C(3)^n2/factorial(n2)*exp(-C(3))*(C(3)+C(1)*s2)^m2/factorial(m2)*exp(-C(3)-C(1)*s2)) ...
%                              ,x0,A,[0,0,0],[],[],lb,[],[],options);
%         MU_f(i,y)=X(1);
%         B1_f(i,y)=X(2);
%         B2_f(i,y)=X(3);
%         Lmin_f(i,y)=fval;
%         [y,N-i]
%     end

    %% Fit mu_hat=f(b_hh)
    mu0=7; options=optimoptions('fminunc','Display','off');
    for i=1:N
        clc
        n1=In1(i,y); m1=Im1(i,y); n2=In2(i,y); m2=Im2(i,y);
        [X,fval]=fminunc(@(P)m1/2+m2/2+n1/2+n2/2 ...
                      -log((1/2^(2*m1)/2^(2*m2)/2^(2*n1)/2^(2*n2)*(m1+n1+2*P*s1+(4*P^2*s1^2-4*P*m1*s1+4*P*n1*s1+m1^2+2*m1*n1+n1^2)^(1/2))^m1) ...
                      /(factorial(m1)*factorial(m2)*factorial(n1)*factorial(n2))) ...
                      -log((m2+n2+2*P*s2+(4*P^2*s2^2-4*P*m2*s2+4*P*n2*s2+m2^2+2*m2*n2+n2^2)^(1/2))^m2) ...
                      -log((m1+n1-2*P*s1+(4*P^2*s1^2-4*P*m1*s1+4*P*n1*s1+m1^2+2*m1*n1+n1^2)^(1/2))^n1) ...
                      -log((m2+n2-2*P*s2+(4*P^2*s2^2-4*P*m2*s2+4*P*n2*s2+m2^2+2*m2*n2+n2^2)^(1/2))^n2) ...
                      +(4*P^2*s1^2-4*P*m1*s1+4*P*n1*s1+m1^2+2*m1*n1+n1^2)^(1/2)/2 ...
                      +(4*P^2*s2^2-4*P*m2*s2+4*P*n2*s2+m2^2+2*m2*n2+n2^2)^(1/2)/2 ...
                    ,mu0,options);

        MU_fm(i,y)=X;
        B1_fm(i,y)=m1/4+n1/4-(X*s1)/2+(4*X^2*s1^2-4*X*m1*s1+4*X*n1*s1+m1^2+2*m1*n1+n1^2)^(1/2)/4;
        B2_fm(i,y)=m2/4+n2/4-(X*s2)/2+(4*X^2*s2^2-4*X*m2*s2+4*X*n2*s2+m2^2+2*m2*n2+n2^2)^(1/2)/4;
        Lmin_fm(i,y)=fval;

        if((n1==0&&m2==0)&&(n2~=0||m1~=0)), comp_ffm(i,y)=-1;
        elseif((m1==0&&n2==0)&&(n1~=0||m2~=0)), comp_ffm(i,y)=1;
        end
%         if((abs(MU_f(i,y)-MU_fm(i,y))<=mu/100)&&(abs(B1_f(i,y)-B1_fm(i,y))<=b1/100)&&(abs(B2_f(i,y)-B2_fm(i,y))<=b2/100)), comp_ffm(i,y)=0;
%         end
    end
    
%     check_ffm(y)=any(isnan(comp_ffm-ones(N,sz)));

    %% Fix for ni=mj=0
%     MUfix=MU_fm; B1fix=B1_fm; B2fix=B2_fm;
%     MUfixUP=MU_fm; B1fixUP=B1_fm; B2fixUP=B2_fm;
%     MUfixDO=MU_fm; B1fixDO=B1_fm; B2fixDO=B2_fm;
%     
%     for i=1:sz
%         for j=1:N
%             if(s1==s2 && In1(j,i)==0 && Im2(j,i)==0)
%                 MUfix(j,i)=(Im1(j,i)-In2(j,i))/4/s1;
%                 B1fix(j,i)=(Im1(j,i)+In2(j,i))/4;
%                 MUfixUP(j,i)=Im1(j,i)/2/s1;
%                 B1fixUP(j,i)=0;
%                 MUfixDO(j,i)=-In2(j,i)/2/s1;
%                 B1fixDO(j,i)=(Im1(j,i)+In2(j,i))/2;
%             elseif(s1==s2 && In2(j,i)==0 && Im1(j,i)==0)
%                 MUfix(j,i)=(Im2(j,i)-In1(j,i))/4/s2;
%                 B2fix(j,i)=(Im2(j,i)+In1(j,i))/4;
%                 MUfixUP(j,i)=Im2(j,i)/2/s2;
%                 B2fixUP(j,i)=0;
%                 MUfixDO(j,i)=-In1(j,i)/2/s2;
%                 B2fixDO(j,i)=(Im2(j,i)+In1(j,i))/2;
%             end
%         end
%     end
% Bfix=(B1fix+B2fix)/2; BfixUP=(B1fixUP+B2fixUP)/2; BfixDO=(B1fixDO+B2fixDO)/2; 
    %% Estimating variance
%     ddLmu=(Im1*s1^2)./(b1+repmat(mu,N,1).*s1).^2+(Im2*s2^2)./(b2+repmat(mu,N,1).*s2).^2;
%     ddLb1=In1/b1^2+Im1./(b1+repmat(mu,N,1).*s1).^2;
%     ddLb2=In2/b2^2+Im2./(b2+repmat(mu,N,1).*s2).^2;
%     ddLmub1=Im1*s1./(b1+repmat(mu,N,1).*s1).^2;
%     ddLmub2=Im2*s2./(b2+repmat(mu,N,1).*s2).^2;
%     dmm=mean(ddLmu);dbb1=mean(ddLb1);dbb2=mean(ddLb2);dmb1=mean(ddLmub1);dmb2=mean(ddLmub2);
%     for i=1:sz
%         I=[dmm(i),dmb1(i),dmb2(i);dmb1(i),dbb1(i),0;dmb2(i),0,dbb2(i)];
%         V=inv(I);
%         vtMU(i)=V(1,1); vtB1(i)=V(2,2); vtB2(i)=V(3,3); vtB(i)=(vtB1(i)+vtB2(i))/2;
%     end

    ddLmu=(Im1*s1^2)./(B1_fm+MU_fm*s1).^2+(Im2*s2^2)./(B2_fm+MU_fm*s2).^2;
    ddLb1=In1/B1_fm^2+Im1./(B1_fm+MU_fm*s1).^2;
    ddLb2=In2/B2_fm^2+Im2./(B2_fm+MU_fm*s2).^2;
    ddLmub1=Im1*s1./(B1_fm+MU_fm*s1).^2;
    ddLmub2=Im2*s2./(B2_fm+MU_fm*s2).^2;
    dmm=mean(ddLmu);dbb1=mean(ddLb1);dbb2=mean(ddLb2);dmb1=mean(ddLmub1);dmb2=mean(ddLmub2);
    if(Im1==0), dmb1=0; dmm=mean((Im2*s2^2)./(B2_fm+MU_fm*s2).^2); end
    if(Im2==0), dmb2=0; dmm=mean((Im1*s1^2)./(B1_fm+MU_fm*s1).^2); end
    I=[dmm,dmb1,dmb2;dmb1,dbb1,0;dmb2,0,dbb2];
    V=inv(I); %covariance matrix
    vtMU=V(1,1); vtB1=V(2,2); vtB2=V(3,3); vtB=(vtB1+vtB2)/2;
    errMU=sqrt(vtMU); errB1=sqrt(vtB1); errB2=sqrt(vtB2);
    
    if(In1==0&&Im1==0)
        I=[dmm,dmb2;dmb2,dbb2]; 
        V=inv(I);
        vtMU=V(1,1); vtB1=Inf; vtB2=V(2,2); vtB=(vtB1+vtB2)/2;
        errMU=sqrt(vtMU); errB1=sqrt(vtB1); errB2=sqrt(vtB2);
    end
    if(In2==0&&Im2==0)
        I=[dmm,dmb1;dmb1,dbb1]; 
        V=inv(I);
        vtMU=V(1,1); vtB1=V(2,2); vtB2=Inf; vtB=(vtB1+vtB2)/2;
        errMU=sqrt(vtMU); errB1=sqrt(vtB1); errB2=sqrt(vtB2);    
    end


    %% Analytical min of CR 
% %     when m1=0 and s1!=2
%     B_1=n1/2;
%     B_2=s2*n2/(s2-s1);
%     MU_=m2/(s1+s2)-n2/(s2-s1);
%     Lmin_=-log(B_1.^n1.*(B_1+MU_*s1).^m1.*exp(-2*B_1-MU_*s1)./(factorial(n1).*factorial(m1))) ...
%            -log(B_2.^n2.*(B_2+MU_*s2).^m2.*exp(-2*B_2-MU_*s2)./(factorial(n2).*factorial(m2)));
% %     when n1=0 and s1!=2
%     B_1=m1/2+s1*(n2/(s1+s2)-m2/(s2-s1));
%     B_2=s2*n2/(s1+s2);
%     MU_=m2/(s2-s1)-n2/(s1+s2);
%     Lmin_=L2bin(B_1,B_2,MU_,s1,s2,n1,m1,n2,m2);
    
end

%%
clear('m1','m2','n1','n2','i','mm1','mm2','nn1','nn2','X','fval','mu0','lb','ub','options','A','x0','ddLmu','ddLb1','ddLb2','ddLmub1','ddLmub2','dmb1','dmb2','dmm','dbb1','dbb2','I','sz','V','vtB','vtB1','vtB2','vtMU','y');
clear('b1','b2','comp_ffm','Im1','Im2','In1','In2','Lmin_fm','mu','N','s1','s2');