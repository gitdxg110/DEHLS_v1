function [GlobalMin, outcome,outeval, CPUtime]= DE_hadamard_CEC2005(fhd,Dimension,Popsize,Maxfeval,VRmin,VRmax,o, A, M, a, alpha, b,varargin)
% 2017-3-27 programming by Xiaogang Dong at jju
    tic;
    global initial_flag
    initial_flag = 0;
    D  = Dimension;
    NP = Popsize;    
    F = 0.9;
    CR = 0.9;
    lowbound  = ones(1,D).*VRmin;
    highbound = ones(1,D).*VRmax;
    pop       = zeros(NP,D);
    outcome=[];
    outeval=[];
    for i=1:NP
        pop(i,:) = lowbound + rand(1,D).*(highbound - lowbound);
    end
    val       = zeros(1,NP);          % create and reset the "cost array"
    nfeval    = 0;                    % number of function evaluations
    for i=1:NP                        % check the remaining members
        val(i)  = feval(fhd,pop(i,:),varargin{:},o, A, M, a, alpha, b);
    end
    [best,bestindex]=min(val);
    nfeval  = nfeval + NP;
    outcome=[outcome best];
    outeval=[outeval nfeval];
    r       = zeros(1,3);  
    iter=2;
    while nfeval <=Maxfeval

        for i=1:NP          
            rd=randperm(NP);
            for j=1:4
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
                    if (tempx(j)<lowbound(j)|| tempx(j)>highbound(j))
                        tempx(j) = lowbound(j) + rand*(highbound(j) - lowbound(j));
                    end
                end
            end
            tempval = feval(fhd, tempx, varargin{:},o, A, M, a, alpha, b);
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
                    for k=1:4
                      if rand<0.5
                        hchild_v=feval(fhd, hchild, varargin{:},o, A, M, a, alpha, b);
                        nfeval  = nfeval + 1;
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
                  
                end
            end            
        end %--for i=1:NP
        outcome=[outcome best];
        outeval=[outeval nfeval];
        iter=iter+1;
    end %--while nfeval <nfevalmax     
    %gbest=bestindex;
    GlobalMin=best;     
    CPUtime   = toc;
end


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
  h=hadamard(8);
  rD=randperm(D);
  r=sort(rD(1:8));   
  s1=[0,r(1),r(2),r(3),r(4),r(5),r(6),r(7),r(8)];
  s2=[r(1),r(2),r(3),r(4),r(5),r(6),r(7),r(8),D];
  for i=1:8   
      for j=1:8
          if h(i,j)==1
            hoffspring4(i,(s1(j)+1):s2(j))=x1((s1(j)+1):s2(j));  
          else
            hoffspring4(i,(s1(j)+1):s2(j))=x2((s1(j)+1):s2(j));   
          end
      end
  end 
end