function [Jsl,Jdl]=cplexDoubleSL(model,cutoff,eliList,atpm)
%% [Jsl,Jdl]=doubleSL(model,cutoff,eliList,atpm)
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
% is true.
% atpm           ATPM Reaction Id in model.rxns if other than 'ATPM'
%OUTPUT
% Jsl        Indices of single lethal reactions identified
% Jdl        Indices of double lethal reactions identified
%
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
        eliList = model.rxns(ismember(model.rxns,'ATPM')); %To eliminate ATPM.
    end
else
    eliList = model.rxns(ismember(model.rxns,'ATPM'));
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

% Wild type solution
solWT=modeldel.solve();
grWT=solWT.objval;
Jnz=find(~eq(solWT.x,0));
if (~isempty(eliList))
    eliIdx = find(ismember(model.rxns,eliList));
    Jnz=Jnz(~ismember(Jnz,eliIdx)); %Jnz
end

Jsl=cplexSingleSL(model,cutoff,eliList);
Jsl=find(ismember(model.rxns,Jsl));

Jnz_copy=Jnz(~ismember(Jnz,Jsl)); %Jnz-Jsl

Jdl=[];

%%
h = waitbar(0,'0.00','Name','Identifying Jdl - Part 1 of 2...');

%Identify Double Lethal Reaction Deletions...
for iRxn=1:length(Jnz_copy)
    delIdx_i=Jnz_copy(iRxn);
    modeldel.Model.lb(delIdx_i)=0;    modeldel.Model.ub(delIdx_i)=0;
    modeldel.Param.simplex.display.Cur=0;
    solKO_i=modeldel.solve();
    newnnz=find(~eq(solKO_i.x,0));
    Jnz_i=newnnz(~ismember(newnnz,Jnz));
    
    if (~isempty(eliList))
        Jnz_i=Jnz_i(~ismember(Jnz_i,eliIdx));
    end
    
    for jRxn=1:length(Jnz_i)
        delIdx_j=Jnz_i(jRxn);
        modeldel.Model.lb(delIdx_j)=0;  modeldel.Model.ub(delIdx_j)=0;
        solKO_ij=modeldel.solve();
        if (solKO_ij.objval<cutoff*grWT || isnan(solKO_ij.objval))
            Jdl=[Jdl;delIdx_i delIdx_j];
        end
        %Reset bounds on idx1 reaction
        modeldel.Model.lb(delIdx_j)=model.lb(delIdx_j);
        modeldel.Model.ub(delIdx_j)=model.ub(delIdx_j);
    end
    %Reset bounds on idx reaction
    modeldel.Model.lb(delIdx_i)=model.lb(delIdx_i);
    modeldel.Model.ub(delIdx_i)=model.ub(delIdx_i);
    waitbar(iRxn/length(Jnz_copy),h,[num2str(round(iRxn*100/length(Jnz_copy))) '% completed...']);

end
close(h);

%%
h = waitbar(0,'0.00','Name','Identifying Jdl - Part 2 of 2...');

for iRxn=1:length(Jnz_copy)
    for jRxn=1:length(Jnz_copy)
        if (jRxn<iRxn)
            delIdx_i=Jnz_copy(iRxn);
            delIdx_j=Jnz_copy(jRxn);
            modeldel.Model.lb(delIdx_i)=0;modeldel.Model.ub(delIdx_i)=0;modeldel.Model.lb(delIdx_j)=0;modeldel.Model.ub(delIdx_j)=0;
            modeldel.Param.simplex.display.Cur=0;
            solKO_ij=modeldel.solve();
            if (solKO_ij.objval<cutoff*grWT || isnan(solKO_ij.objval))
                Jdl=[Jdl;delIdx_i delIdx_j];
            end
            modeldel.Model.lb(delIdx_i)=model.lb(delIdx_i);  modeldel.Model.ub(delIdx_i)=model.ub(delIdx_i);
            modeldel.Model.lb(delIdx_j)=model.lb(delIdx_j);  modeldel.Model.ub(delIdx_j)=model.ub(delIdx_j);
        else
            break;
        end
    end
    waitbar(iRxn*(iRxn-1)/(length(Jnz_copy)*(length(Jnz_copy)-1)),h,[num2str(round(iRxn*(iRxn-1)*100/(length(Jnz_copy)*(length(Jnz_copy)-1)))) '% completed...']);

end
close(h);
Jsl=model.rxns(Jsl);
Jdl=model.rxns(Jdl);

fprintf('\n Done...');
end