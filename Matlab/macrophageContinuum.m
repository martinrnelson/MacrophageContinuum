clear

%% Timestepping parameters

tMax=1000;
tSpan=[0 tMax];

%% Discretisation of phenotypes (p)

nPhenoPoints=101;
p=linspace(-1,1,nPhenoPoints);
p=p';
dp=p(2)-p(1);

%% Model Parameters

% Set parameters to the values of Table 1
params=setBaselineParams(); 

% Place modifications to the standard parameter set below, if needed
params.gammaG=0.2;          % To reproduce Figure 2(b), for example

%% Initial conditions

g0=0.5;
c0=0.5;
m0=0.*p+10;

ICs=[g0;c0;m0];

%% Mediator production functions

f1=@(p)((1-p)/2);
f2=@(p)((p+1)/2);  

%% Macrophage recruitment function (Choose either R=1 or R=Gaussian)

% For R=1
R=@(p)(0.*p+1);   

% For R=Gaussian (mean mu, s.d. sigma)
% G=@(p,mu,sigma)(exp(-(p-mu).^2/sigma^2));
% R=@(p)(G(p,-1,0.1));

%% Create the two finite difference matrices (A), using upwinding
Afwd=zeros(nPhenoPoints);

for i=1:nPhenoPoints
   Afwd(i,i)=-1;
end
for i=1:nPhenoPoints-1
   Afwd(i,i+1)=1;
end
Afwd=Afwd/dp;

Aback=zeros(nPhenoPoints);
for i=1:nPhenoPoints
   Aback(i,i)=1;
end
for i=2:nPhenoPoints
  Aback(i,i-1)=-1;
end  
Aback=Aback/dp;

%% ODE Solver

% Solver options
opts=odeset('RelTol',1e-10,'AbsTol',1e-12);

% Call ODE solver
sol = ode45(@(t,y)macrophageContinuum_RHS(t,y,p,Afwd,Aback,f1,f2,R,params),tSpan,ICs,opts);

% Store results
t=sol.x;
g=sol.y(1,:);
c=sol.y(2,:);
m=sol.y(3:end,:);

% Compute total macrophage number (as function of t)
for j=1:numel(t)
    M(j)=trapz(p,m(:,j));    
end

% Compute the median phenotype (as a function of t)
pHat=zeros(size(t));
cumulativeMacros=cumsum(m); % This sums across phenotypes for fixed t
for j=1:numel(t)
    idx=find(cumulativeMacros(:,j)>M(j)/2,1,'first');
    if(numel(idx)<1)
        pHat(j)=-2;
    else
        pHat(j)=p(idx);
    end
end

%% Plot solutions

close all
figure(1);
plot(t,g,'b');
hold on
plot(t,c,'r');
plot(t,M,'g');
xlabel('t');

figure(2);
imagesc(t,p,m);
set(gca,'YDir','normal');
xlabel('t');
ylabel('p');
caxis([0 25]);
colorbar;

figure(3);
plot(t,pHat,'k');
ylim([-2 2]);