function [gbest,gbestval,fitcount]= DE_hadamard_TEC13(fhd,Dimension,Popsize,maxFES,VRmin,VRmax,time,varargin)
% 2014-3-25 programming by Hu Peng at WHU
    tic;
   
    global initial_flag
    initial_flag = 0;
    D  = Dimension;
    NP = Popsize;    
    F = 0.9;
    CR = 0.9; 
    trace=[];
    t_FES=[];
    lowbound  = ones(1,D).*VRmin;
    highbound = ones(1,D).*VRmax;
    pop       = zeros(NP,D);
    for i=1:NP
        pop(i,:) = lowbound + rand(1,D).*(highbound - lowbound);
    end
    val       = zeros(1,NP);          % create and reset the "cost array"
    nfeval    = 0;                    % number of function evaluations
    for i=1:NP                        % check the remaining members
        val(i)  = feval(fhd,pop(i,:),varargin{:});
    end
    [best,bestindex]=min(val);
    nfeval  = nfeval + NP;
    r       = zeros(1,3);
    trace=[trace best];
    t_FES=[t_FES nfeval];
    iter=2;
    while nfeval <=maxFES

        for i=1:NP
            rd=randperm(NP);
            for j=1:3
                if rd(j)~=i 
                    r(j)=rd(j);
                else
                    r(j)=rd(5);
                end
            end 
            tempx = pop(i, :);            
            jrand=fix(1+D*rand);
            for j=1:D
                if (rand<CR || j==jrand ) 
                    tempx(j) = pop(r(1),j)+F.*(pop(r(2),j)-pop(r(3),j));
                    %tempx(j) = pop(bestindex,j)+F.*(pop(r(1),j)-pop(r(2),j));
                    %tempx(j) = pop(r(1),j)+F.*(pop(bestindex,j)-pop(r(1),j))+F.*(pop(r(2),j)-pop(r(3),j));
                    %tempx(j) = pop(i,j)+F.*(pop(bestindex,j)-pop(i,j))+F.*(pop(r(1),j)-pop(r(2),j));
                    if (tempx(j)<lowbound(j)|| tempx(j)>highbound(j))
                        tempx(j) = lowbound(j) + rand*(highbound(j) - lowbound(j));
                    end
                end
            end
            tempval = feval(fhd, tempx, varargin{:});
            nfeval  = nfeval + 1;
            if tempval < val(i)
                val(i)   = tempval;
                pop(i,:) = tempx;
                 if tempval<=best
                      best=tempval;
                      bestindex=i;
                 end
            else
                if rand<0.1
                    hchild=hadamard_func7(tempx,pop(i,:),D);
                    hchild_v=feval(fhd, hchild, varargin{:});
                    nfeval  = nfeval + 4;
                    [hbest,hbest_index]=min(hchild_v);
                    if hbest < val(i)
                         val(i)   = hbest;
                         pop(i,:) = hchild(hbest_index);
                        if hbest<=best
                           best=hbest;
                           bestindex=i;
                        end
                    end
                end
            end            
        end %--for i=1:NP
        trace=[trace best];
        t_FES=[t_FES nfeval];
        iter=iter+1;      

    end %--while nfeval <nfevalmax
    filename=['E:\Rearch_Program\hadamardDE1\5new_DE_hadamard30P\TEC\trace\trace_',num2str(varargin{:}),'_',num2str(time),'.mat'];
    save(filename,'trace','t_FES');   
    gbest=bestindex;
    gbestval=best; 
    fitcount=nfeval;
    CPUtime   = toc;
end

% 
% function hoffspring3=hadamard_func(x1,x2,D)
%   h=hadamard(4);
%   rD=randperm(D);
%   r=sort(rD(1:3));  
%   for i=2:4   
%       if h(i,1)==1
%           hoffspring3(i-1,1:r(1))=x1(1:r(1));
%       else
%           hoffspring3(i-1,1:r(1))=x2(1:r(1));
%       end
%       if h(i,2)==1
%           hoffspring3(i-1,(r(1)+1):r(2))=x1((r(1)+1):r(2));
%       else
%           hoffspring3(i-1,(r(1)+1):r(2))=x2((r(1)+1):r(2));
%       end
%       if h(i,3)==1
%           hoffspring3(i-1,(r(2)+1):r(3))=x1((r(2)+1):r(3));
%       else
%           hoffspring3(i-1,(r(2)+1):r(3))=x2((r(2)+1):r(3));
%       end
%       if h(i,4)==1
%         hoffspring3(i-1,(r(3)+1):D)=x1((r(3)+1):D);
%       else
%         hoffspring3(i-1,(r(3)+1):D)=x2((r(3)+1):D); 
%       end
%   end
% end
function hoffspring4=hadamard_func7(x1,x2,D)
  h=hadamard(4);
  rD=randperm(D);
  r=sort(rD(1:3));   
  s1=[0,r(1),r(2),r(3)];
  s2=[r(1),r(2),r(3),D];
  for i=1:4   
      for j=1:4
          if h(i,j)==1
            hoffspring4(i,(s1(j)+1):s2(j))=x1((s1(j)+1):s2(j));  
          else
            hoffspring4(i,(s1(j)+1):s2(j))=x2((s1(j)+1):s2(j));   
          end
      end
  end 
end