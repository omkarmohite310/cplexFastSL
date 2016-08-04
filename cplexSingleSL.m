function [Jsl]=cplexSingleSL(model,cutoff,eliList,atpm)
%% [Jsl]=singleSL(model,cutoff,eliList,atpm)
% INPUT
% model (the following fields are required - others can be supplied)       
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%   rxns         Reaction Names
%OPTIONAL
% cutoff         cutoff percentage value for lethality.Default is 0.01.
% eliList        List of reactions to be ignored for lethality
% analysis:Exchange Reactions, ATPM etc.
% atpm           ATPM Reaction Id in model.rxns if other than 'ATPM'
%OUTPUT
% Jsl            Single lethal reactions identified
% Aditya Pratapa       6/26/14. 

if exist('cutoff', 'var')
    if isempty(cutoff)
        cutoff = 0.01;
    end
else
    cutoff = 0.01;
end

if exist('atpm', 'var')
    if isempty(atpm)
        atpm = 'ATPM'; %Reaction Id of ATP maintenance reaction- by default it takes 'ATPM'
    end
else
    atpm = 'iAF_ATPM';
end



if exist('eliList', 'var')
    if isempty(eliList)
        eliList = model.rxns(ismember(model.rxns,atpm)); %To eliminate ATPM.
    end
else
        eliList = model.rxns(ismember(model.rxns,atpm));
end

% Defining CPlex optimization problem
[nMets,nRxns]=size(model.S);
modeldel=Cplex();
modeldel.Model.A=sparse(model.S);
modeldel.Model.obj=model.c;
modeldel.Model.rhs=model.b;
modeldel.Model.lhs=model.b;
modeldel.Model.lb=model.lb;
modeldel.Model.ub=model.ub;
modeldel.Model.sense='maximize';
modeldel.Param.barrier.display.Cur=0;
modeldel.Param.simplex.display.Cur=0;

Jsl=[];

%Step1 Identify Single Lethal Reactions...
%Identify minNorm flux distribution
solWT=modeldel.solve();
grWT=solWT.objval;

Jnz=find(~eq(solWT.x,0));

if (~isempty(eliList))
    eliIdx = find(ismember(model.rxns,eliList));
    Jnz=Jnz(~ismember(Jnz,eliIdx)); %Jnz
end
h = waitbar(0,'0.00','Name','Identifying Jsl...');

%Identify Single Lethal Reaction Deletions...
for iRxn=1:length(Jnz)
    delIdx_i=Jnz(iRxn);
    modeldel.Model.lb(delIdx_i)=0;    modeldel.Model.ub(delIdx_i)=0;
    modeldel.Param.simplex.display.Cur=0;
    solKO_i=modeldel.solve();
    if (solKO_i.objval<cutoff*grWT || isnan(solKO_i.objval))
        Jsl=[Jsl;delIdx_i];
    end
   %Reset bounds on idx reaction
   modeldel.Model.lb(delIdx_i)=model.lb(delIdx_i);
   modeldel.Model.ub(delIdx_i)=model.ub(delIdx_i);
   waitbar(iRxn/length(Jnz),h,[num2str(round(iRxn*100/length(Jnz))) '% completed...']);
    
end
close(h);
Jsl=model.rxns(Jsl);
end