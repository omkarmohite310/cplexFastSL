function [Jsl,Jdl,Jtl]=cplexTripleSL(model,cutoff,eliList,atpm)
%%  [slist_id,dlist_id,tlist_id]=tripleSL(model,cutoff,eliList,atpm)
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
% Jtl        Indices of triple lethal reactions identified
% Aditya Pratapa       7/1/14.
%%

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
    atpm = 'ATPM';
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


%Wildtype solution
solWT=modeldel.solve();
grWT=solWT.objval;
Jnz=find(~eq(solWT.x,0));
%If a list of reactions for which are eliminated for lethality is given often exchange reactions are not considered
if (~isempty(eliList))
    eliIdx = find(ismember(model.rxns,eliList));
    Jnz_copy=Jnz(~ismember(Jnz,eliIdx)); %Jnz
end

Jsl=cplexSingleSL(model,cutoff,eliList);
Jsl=find(ismember(model.rxns,Jsl));
Jnz_copy=Jnz_copy(~ismember(Jnz_copy,Jsl)); %Eliminate Single lethal reaction deletions for enumeration of higher order lethals

Jdl=[];
Jtl=[];
%%

h = waitbar(0,'0.00','Name','Identifying Jdl & Jtl - Part 1 of 2...');

%Identify Double and Triple Lethal Reaction Deletions...
for iRxn=1:length(Jnz_copy)
    delIdx_i=Jnz_copy(iRxn);
    modeldel.Model.lb(delIdx_i)=0;    modeldel.Model.ub(delIdx_i)=0;
    modeldel.Param.simplex.display.Cur=0;
    solKO_i=modeldel.solve();  %It can't be a single lethal so we can proceed further
    newnnz=find(~eq(solKO_i.x,0));
    Jnz_i=newnnz(~ismember(newnnz,Jnz));
    
    if (~isempty(eliList))
        Jnz_i=Jnz_i(~ismember(Jnz_i,eliIdx)); %Eliminate Exchange and ATP Maintenance reactions
    end
    
    for jRxn=1:length(Jnz_i)
            delIdx_j=Jnz_i(jRxn);
            modeldel.Model.lb(delIdx_j)=0;modeldel.Model.ub(delIdx_j)=0;
            solKO_ij=modeldel.solve();
            if (solKO_ij.objval<cutoff*grWT && ~eq(solKO_ij.status,0)) 
                Jdl=[Jdl;delIdx_i delIdx_j];
            else
                  if eq(solKO_ij.status,0)
                    solKO_ij=modeldel.solve();
                    if (solKO_ij.objval<cutoff*grWT || isnan(solKO_ij.objval)) 
                        Jdl=[Jdl;delIdx_i delIdx_j];
                        modeldel.Model.lb(delIdx_j)=model.lb(delIdx_j);
                        modeldel.Model.ub(delIdx_j)=model.ub(delIdx_j);
                        continue;
                    end
                end
                Jnz_ij=find(~eq(solKO_ij.x,0));
                Jnz_ij=Jnz_ij(~ismember(Jnz_ij,Jnz));
                
                if (~isempty(eliList))
                    Jnz_ij=Jnz_ij(~ismember(Jnz_ij,eliIdx)); %Eliminate Exchange and ATPM reactions
                end
                
                for kRxn=1:length(Jnz_ij)           
                    
                    delIdx_k=Jnz_ij(kRxn);
                    modeldel.Model.lb(delIdx_k)=0;modeldel.Model.ub(delIdx_k)=0;
                    solKO_ijk=modeldel.solve();
                    if (solKO_ijk.objval<cutoff*grWT || isnan(solKO_ijk.objval))
                        Jtl=[Jtl;delIdx_i delIdx_j delIdx_k];
                    end
                    
                    modeldel.Model.lb(delIdx_k)=model.lb(delIdx_k);
                    modeldel.Model.ub(delIdx_k)=model.ub(delIdx_k);             
                end
            end
            modeldel.Model.lb(delIdx_j)=model.lb(delIdx_j);
            modeldel.Model.ub(delIdx_j)=model.ub(delIdx_j);
       
    end
    
    modeldel.Model.lb(delIdx_i)=model.lb(delIdx_i);
    modeldel.Model.ub(delIdx_i)=model.ub(delIdx_i);
    waitbar(iRxn/length(Jnz_copy),h,[num2str(round(iRxn*100/length(Jnz_copy))) '% completed...']);

end
close(h);

%%
h = waitbar(0,'0.00','Name','Identifying Jdl & Jtl - Part 2 of 2...');

for iRxn=1:length(Jnz_copy)
    for jRxn=1:length(Jnz_copy)
        if (jRxn<iRxn)
            delIdx_i=Jnz_copy(iRxn);
            delIdx_j=Jnz_copy(jRxn);
            modeldel.Model.lb(delIdx_i)=0;modeldel.Model.ub(delIdx_i)=0;
            modeldel.Model.lb(delIdx_j)=0;modeldel.Model.ub(delIdx_j)=0;
            solKO_ij=modeldel.solve();
            if (solKO_ij.objval<cutoff*grWT && ~eq(solKO_ij.status,0))
                Jdl=[Jdl;delIdx_i delIdx_j];
            else
                  if eq(solKO_ij.status,0)
                    solKO_ij=modeldel.solve();
                    if (solKO_ij.objval<cutoff*grWT || isnan(solKO_ij.objval)) 
                        Jdl=[Jdl;delIdx_i delIdx_j];
                        modeldel.Model.lb(delIdx_j)=model.lb(delIdx_j);
                        modeldel.Model.ub(delIdx_j)=model.ub(delIdx_j);
                        continue;
                    end
                  end
                Jnz_ij=find(~eq(solKO_ij.x,0));
                Jnz_ij=Jnz_ij(~ismember(Jnz_ij,Jnz));
                
                if (~isempty(eliList))
                    Jnz_ij=Jnz_ij(~ismember(Jnz_ij,eliIdx)); %Eliminate Exchange and ATPM reactions
                end
                
                for kRxn=1:length(Jnz_ij)
                    delIdx_k=Jnz_ij(kRxn);
                    modeldel.Model.lb(delIdx_k)=0;modeldel.Model.ub(delIdx_k)=0;
                    solKO_ijk=modeldel.solve();
                    if (solKO_ijk.objval<cutoff*grWT || isnan(solKO_ijk.objval))
                        Jtl=[Jtl;delIdx_i delIdx_j delIdx_k];
                    end
                    
                    modeldel.Model.lb(delIdx_k)=model.lb(delIdx_k);
                    modeldel.Model.ub(delIdx_k)=model.ub(delIdx_k);
                    
                end
                
                for kRxn=1:length(Jnz_copy)
                    
                    if (kRxn<jRxn)
                        delIdx_k=Jnz_copy(kRxn);
                        modeldel.Model.lb(delIdx_k)=0;modeldel.Model.ub(delIdx_k)=0;
                        solKO_ijk=modeldel.solve();
                        if (solKO_ijk.objval<cutoff*grWT ||isnan(solKO_ijk.objval))
                            Jtl=[Jtl;delIdx_i delIdx_j delIdx_k ];
                        end
                        
                        modeldel.Model.lb(delIdx_k)=model.lb(delIdx_k);
                        modeldel.Model.ub(delIdx_k)=model.ub(delIdx_k);
                        
                    else
                        break;
                    end
                end
            end
            modeldel.Model.lb(delIdx_j)=model.lb(delIdx_j);
            modeldel.Model.ub(delIdx_j)=model.ub(delIdx_j);
        else
            break;
        end
        
    end
    modeldel.Model.lb(delIdx_i)=model.lb(delIdx_i);
    modeldel.Model.ub(delIdx_i)=model.ub(delIdx_i);
    waitbar(iRxn*(iRxn-1)*(iRxn-2)/(length(Jnz_copy)*(length(Jnz_copy)-1)*(length(Jnz_copy)-2)),h,[num2str(round(iRxn*(iRxn-1)*(iRxn-2)*100/(length(Jnz_copy)*(length(Jnz_copy)-1)*(length(Jnz_copy)-2)))) '% completed...']);

end
close(h);

%Eliminate double lethal reaction deletions in triple lethal reactions
temporary=[];
g=zeros(1,length(Jdl));
for iRxn=1:length(Jtl)
    for jRxn=1:length(Jdl)
        g(jRxn)=sum(ismember(Jtl(iRxn,:),Jdl(jRxn,:)));
        if g(jRxn)>=2
            break;
        end
    end
    if max(g)<2
        temporary=[temporary;Jtl(iRxn,:)];
    end
end
Jtl=temporary;

%Eliminate duplicates in triple reaction deletions
Jtl=unique(sort(Jtl,2),'rows');



Jsl=model.rxns(Jsl);
Jdl=model.rxns(Jdl);
Jtl=model.rxns(Jtl);
end

