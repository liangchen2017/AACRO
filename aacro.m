

function [fxmin,xmin,Swarm,history] = cro(croOptions,Loninput)
 global fSwarm ;%x1 x2 x3 x4 x5 x6 x7 x8 x9 x10
 global MinHit;
 global d;
 global h;
 global f;
 global buffer;
 global popsize;
 global MinStruct;
 global StepSize;
 global NumHit;
 global PBest;
 global VStep;
 global w_now;
 global Vmax;
 global decom;
 global onwall;
 global inter;
 global synthe;
 decom=0;
 onwall=0;
 inter=0;
 synthe=0;
 buffer=0;

upbnd = 300; % Upper bound for init. of the swarm
lwbnd = 0; % Lower bound for init. of the swarm
GM = 0; % Global minimum (used in the stopping criterion)
ErrGoal = 1e-1; % Desired accuracy
% 

%Initializations
% w_start = croOptions.PSO.w_start;   %Initial inertia weight's value    
% w_end = croOptions.PSO.w_end;       %Final inertia weight    
% w_varyfor = floor(croOptions.PSO.w_varyfor*croOptions.Vars.Iterations);  %Weight change step. Defines total number of iterations for which weight is changed.
% w_now = w_start;
% inertdec = (w_start-w_end)/w_varyfor; %Inertia weight's change per iteration  
    


success = 0; % Success Flag
iter = 0;   % Iterations' counter
fevals = 0; % Function evaluations' counter  ??????????????
%MinStruct=[0 0 0 0 0 0 0 0 0];    

% Initialize Swarm and Velocity
popsize = croOptions.Vars.PopSize;
VirtualSwarm = rand(2*popsize, croOptions.Vars.Dim)*diag(croOptions.Obj.ub-croOptions.Obj.lb) + repmat(croOptions.Obj.lb,2*popsize,1);
En=repmat([0 0],popsize,1);
NumHit=1;
MinHit=1;
StepSize=(croOptions.Obj.ub-croOptions.Obj.lb)/2;
VStep = rand(1, croOptions.Vars.Dim)*diag(croOptions.PSO.Vmax);
Vmax=croOptions.PSO.Vmax;
%1/5 success rule

mutations=croOptions.Vars.Iterations/100;%1/5 success rule eachiternumber
successiters=0;
number=1;%1/5success rule


f2eval = croOptions.Obj.f2eval; %The objective function to optimize.

 fSwarm= feval(f2eval,VirtualSwarm,Loninput); 
 VirtualPopSize=popsize*2;
 for i=1:popsize
     [max1,d]=max(fSwarm);
     if d==VirtualPopSize
         fSwarm(d,:)=[];
         VirtualSwarm(d,:)=[];
     else
         fSwarm(d,:)=fSwarm(VirtualPopSize,:);
         VirtualSwarm(d,:)=VirtualSwarm(VirtualPopSize,:);
         fSwarm(VirtualPopSize,:)=[];
         VirtualSwarm(VirtualPopSize,:)=[];
     end
     VirtualPopSize=VirtualPopSize-1;
 end
 
     Swarm=VirtualSwarm;
     En(:,1)=fSwarm(:,1);
 
%  x1=Swarm(:,1);x2 = Swarm(:,2); x3 = Swarm(:,3); x4 = Swarm(:,4); x5 =Swarm(:,5); x6 = Swarm(:,6); x7 = Swarm(:,7); x8 = Swarm(:,8); x9 =Swarm(:,9); x10 = Swarm(:,10)
[maxpe,x]=max(fSwarm(:,1));
[minpe,y]=min(fSwarm(:,1));
% En(:,2)=maxpe-minpe;

En(:,2)=1000;

% Initializing the MinStruct positions matrix and  
% the corresponding function values 
PBest = Swarm;             
fPBest = fSwarm;          

% Finding MinStruct particle in initial population  
[fGBest, g] = min(fSwarm);
MinPE= fGBest;
MinPE1=MinPE;
MinStruct = Swarm(g,:); %Used to keep track of the MinStruct particle ever
history = [ MinHit,MinPE];

if croOptions.Disp.Interval && (rem(iter, croOptions.Disp.Interval) == 0)
    disp(sprintf('Iterations\t\tfGBest\t\t\tfevals'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  THE  CRO  LOOP                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 tic;
while( (iter < croOptions.Vars.Iterations) )
    iter = iter+1;
  
    % Update the value of the inertia weight w 
%     if ((iter<=w_varyfor) && (iter > 1))
%         w_now = w_now - inertdec; %Change inertia weight 
%     end
    if iter>100000
        a=1;
    end
    MinHitPre=MinHit;
    NumHitpre=NumHit;
    numberpre=number;
if MinPE<50
    a=1;
end
if iter>15000
    a=1;
end
    if(iter>=mutations*number)
        number=number+1;
        successiters(number,1)=0;
    end
    if((number<10)&&(numberpre<number))
       StepSize=StepSize*0.5; 
    end
if ((number>=10)&&(numberpre<number))
    totalsuccess=0;
    for i=1:10
        totalsuccess=totalsuccess+successiters(number-i+1,1);
    end
    if(totalsuccess>2*mutations)
        StepSize=StepSize/0.85;
        if(StepSize>(croOptions.Obj.ub-croOptions.Obj.lb)/2)
            StepSize=(croOptions.Obj.ub-croOptions.Obj.lb)/2;
        end
    else
        if (MinPE>80)
            StepSize=StepSize*0.45;
        else
            StepSize=StepSize*0.5;
        end
    end
        
        
end

    changerate=croOptions.SParams.ChangeRate;
     b=rand(1);
 if b>croOptions.SParams.MoleColl
     c=1+(popsize-1).*rand(1);
     d=round(c);
     if(d>=popsize)
         d=popsize;
     end
     R1=Swarm(d,:); 

     if (iter/croOptions.Vars.Iterations)==0.9
        [MinPE1,Swarm,En]=decomposition(croOptions,f2eval,En,Swarm,MinPE1,Loninput);
     else
        [MinPE1,En,Swarm]=onwallineffectivecollision(croOptions,f2eval,R1,En,Swarm,MinPE1,Loninput);
     end
 else
     e=1+(popsize-1).*rand(1);
     f=round(e);
     if(f>=popsize)
         f=popsize;
     end
     R2=Swarm(f,:);
     g=1+(popsize-1).*rand(1);
     h=round(g);
     if(h>=popsize)
        h=popsize;
     end
     if(f==h)
       if(f==1)
           g=2+(popsize-2).*rand(1);
           h=round(g);
       elseif(f==popsize)
           g=1+(popsize-2).*rand(1);
           h=round(g);
           if(h==0)
            h=1;
           end
       else
        a=rand(1);
        if(a>0.5)
            h=popsize;
        else
            h=1;
        end
       end
      end
     R3=Swarm(h,:);
     temp=rand(1);
     tagpopsize=(croOptions.Vars.PopSize-popsize)/croOptions.Vars.PopSize;
     if((changerate>=temp)&&(tagpopsize<=0.5))
         tag=1;
     else
         tag=0;
     end
     if tag==1
         [MinPE1,En,Swarm]=synthesis(croOptions,f2eval,R2,R3,En,Swarm,MinPE1,Loninput);
     else
         [MinPE1,En,Swarm]=interineffectivecollision(croOptions,f2eval,R2,R3,En,Swarm,MinPE1,Loninput);
     end
 end  
         
   
    if(NumHitpre~=NumHit)  
        successiters(number,1)=successiters(number,1)+1;
    end   
    if iter==croOptions.Vars.Iterations
        [MinPE1,En,Swarm]=final(croOptions,f2eval,Swarm,En,MinPE1,Loninput);
    end
    MinPE=MinPE1;
    %%OUTPUT%%
    if croOptions.Save.Interval && (rem(iter, croOptions.Save.Interval) == 0)&&(MinHitPre~=MinHit)
        history((size(history,1)+1), :) = [MinHit,MinPE];         %SIZE(X,1) returns the number of rows
    end
    
    if croOptions.Disp.Interval && (rem(iter, croOptions.Disp.Interval) == 0)
        disp(sprintf('%4d\t\t\t%.5g\t\t\t%5d\t\t\t%5d\t\t\t%.5g\t\t\t%5d\t\t\t%5d\t\t\t%5d\t\t\t%5d', iter, MinPE, MinHit,popsize,buffer,onwall,decom,synthe,inter));
    end

    
    %%TERMINATION%%
           %zm
    if abs(MinPE) <= croOptions.Vars.ErrGoal     %GBest
        success = 1;
    end

    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  END  OF CRO  LOOP                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 toc


fxmin = MinPE;   %fxmin 
xmin  = MinStruct;
end
 function[MinPE1,Swarm,En]=decomposition(croOptions,f2eval,En,Swarm,MinPE1,Loninput)

    global MinHit;
    global MinStruct;
    global buffer;
    global popsize;
    global StepSize;
    global NumHit;
    global PBest;
    global decom;
    decom=decom+1;
    R4=0;
    R5=0;
    temppopsize=popsize;
    for i=1:temppopsize
        tempStepSize=StepSize;
        maxStepSize=(croOptions.Obj.ub-croOptions.Obj.lb)/2;
        number=1000;
        R1=Swarm(i,:); 
        k=1;
        while(tempStepSize(1,1)<maxStepSize(1,1))
            while number>0
                for j=1:croOptions.Vars.Dim
                    c=tempStepSize(:,j);
                    c1=normrnd(0,c);        
                    R4(:,j)=R1(:,j)+c1;
                    if R4(:,j)<croOptions.Obj.lb(:,j)
                        R4(:,j) = 2*croOptions.Obj.lb(:,j)-R4(:,j);
                    elseif R4(:,j)>croOptions.Obj.ub(:,j)
                        R4(:,j) = 2*croOptions.Obj.ub(:,j)-R4(:,j);
                    end
                    c2=normrnd(0,c); 
                    R5(:,j)=R1(:,j)+c2;
                    if R5(:,j)<croOptions.Obj.lb(:,j)
                        R5(:,j) = 2*croOptions.Obj.lb(:,j)-R5(:,j);
                    elseif R5(:,j)>croOptions.Obj.ub(:,j)
                        R5(:,j) = 2*croOptions.Obj.ub(:,j)-R5(:,j);
                    end
                    e1=rand(1);e2=rand(1);
                    if (e1>=0.6)&&(e2>=0.6)
                        R4(:,j)=R1(:,j);
                    elseif(e1<=0.4)&&(e2<=0.4)
                        R5(:,j)=R1(:,j);
                    end   
                end
                fSwarm1= feval(f2eval,R4,Loninput);
                fSwarm2= feval(f2eval,R5,Loninput);
                A=[fSwarm1;fSwarm2];
                B=[R4;R5];
                [minpe1,x]=min(A);
                [maxpe1,y]=max(A);
                R=B(x,:);
                R1=B(y,:);
                tempmax(k,:)=R1;
                temppemax(k,:)=maxpe1;
                if (minpe1<MinPE1)
                    if(popsize<(croOptions.Vars.PopSize)*1.5)
                        popsize=popsize+1;
                        MinStruct=R;
                        Swarm(popsize,:)=R;
                        En(popsize,1)=minpe1;
                        En(popsize,2)=minpe1;
                        buffer=buffer-minpe1;
                        PBest(popsize,:)=R;
                    else
                        Swarm(i,:)=R;
                        MinStruct=R;
                        buffer=buffer+En(i,1)-minpe1;
                        En(i,1)=minpe1;
                        PBest(i,:)=R;
                    end
                end
                number=number-1;
                k=k+1;
            end 
            tempStepSize=tempStepSize/0.6; 
        end
        [minpe2,xx]=min(temppemax);
        if minpe2<En(i,1)
            Swarm(i,:)=tempmax(xx,:);
            PBest(i,:)=Swarm(i,:);
            buffer=buffer+En(i,1)-minpe2;
            En(i,1)=minpe2;
        end
    end
    NumHit=NumHit+1;
    if popsize>temppopsize
        MinHit=MinHit+1;
    end
 end
 function[MinPE1,En,Swarm]=onwallineffectivecollision(croOptions,f2eval,R1,En,Swarm,MinPE1,Loninput)
   % global  fSwarm ; %x1 x2 x3 x4 x5 x6 x7 x8 x9 x10
    global MinHit;
    global MinStruct;
    global d;
    global buffer;
    global StepSize;
    global NumHit;
    global PBest;
    global onwall;
    onwall=onwall+1;

    R6=0;
    energytemp=En(d,1)+En(d,2);
    r=rand(1);
    if r>=croOptions.SParams.w_global
        tag=1;
    else
        tag=0;
    end

    for i=1:croOptions.Vars.Dim
        c=StepSize(:,i);
        c1=normrnd(0,c);
        if(c1<=-c)
            c1= -c;
        elseif(c1>=c)
            c1=c;
        end
        if tag==1
            r1=rand(1);
            r2=rand(1);
            R6(:,i)=R1(:,i)+c1+croOptions.PSO.c1*r1*(PBest(d,i)-R1(:,i))+croOptions.PSO.c2*r2*(MinStruct(:,i)-R1(:,i));
        else
            R6(:,i)=R1(:,i)+c1; 
        end
        if R6(:,i)<croOptions.Obj.lb(:,i)
            R6(:,i) = 2*croOptions.Obj.lb(:,i)-R6(:,i);
        elseif R6(:,i)>croOptions.Obj.ub(:,i)
        R6(:,i) = 2*croOptions.Obj.ub(:,i)-R6(:,i);
        end
    end
    
     fSwarm3= feval(f2eval,R6,Loninput);
     pe=fSwarm3;
 if(energytemp>pe)
     energy=energytemp-pe;
     LossRate=croOptions.SParams.KELossRate;
     r=LossRate+(1-LossRate).*rand(1);
     ke=energy*r;
     buffer=buffer+energy*(1-r);
     Swarm(d,:)=R6;
     En(d,1)=pe;
     En(d,2)=ke;
     PBest(d,:)=R6;
     NumHit=NumHit+1;
 end
 if(pe<MinPE1)
        MinPE1=pe;
        MinStruct=R6;
        MinHit=NumHit; 
 end


 end
 function[MinPE1,En,Swarm]=synthesis(croOptions,f2eval,R2,R3,En,Swarm,MinPE1,Loninput)
    %global  fSwarm ; %x1 x2 x3 x4 x5 x6 x7 x8 x9 x10
    global MinStruct;
    global MinHit;
    global f;
    global h;
    global popsize;
    global NumHit;
    global PBest;
    global StepSize;
    global synthe;
    synthe=synthe+1;
    R7=0;
    R10=0;
    energytemp=En(f,1)+En(f,2)+En(h,1)+En(h,2);
        for i=1:croOptions.Vars.Dim 
            k=rand(1);
            if(k>=0.5)
             z=1;x=0;
            else
            z=0;x=1; 
            end
            R7(1,i) = z*R2(1,i)+ x*R3(1,i);
            R10(1,i)=x*R2(1,i)+z*R3(1,i);
            c=StepSize(:,i);
            C1(1,:)=normrnd(0,c);
            if(C1(1,:)<=croOptions.Obj.lb(:,i))
                C1(1,:)= croOptions.Obj.lb(:,i);
            elseif(C1(1,:)>=croOptions.Obj.ub(:,i))
                C1(1,:)=croOptions.Obj.ub(:,i);
            end
            C2(1,:)=normrnd(0,c);
            if(C2(1,:)<=croOptions.Obj.lb(:,i))
                C2(1,:)= croOptions.Obj.lb(:,i);
            elseif(C2(1,:)>=croOptions.Obj.ub(:,i))
                C2(1,:)=croOptions.Obj.ub(:,i);
            end
        end
        Q = rand(1, croOptions.Vars.Dim);
        W = rand(1, croOptions.Vars.Dim);
        
        TEMP1 = R7+C1+croOptions.PSO.c1*Q.*(PBest(f,:)-R7) + croOptions.PSO.c2*W.*(MinStruct-R7);
        TEMP2 = R10+C2+croOptions.PSO.c1*Q.*(PBest(h,:)-R10) + croOptions.PSO.c2*W.*(MinStruct-R10);
        changeRow1=TEMP1>croOptions.Obj.ub;
        TEMP1(find(changeRow1))=2*croOptions.Obj.ub(find(changeRow1))-TEMP1(find(changeRow1));
        changeRow2=TEMP1<croOptions.Obj.lb;
        TEMP1(find(changeRow2))=2*croOptions.Obj.lb(find(changeRow2))-TEMP1(find(changeRow2)); 
        changeRow3=TEMP2>croOptions.Obj.ub;
        TEMP2(find(changeRow3))=2*croOptions.Obj.ub(find(changeRow3))-TEMP2(find(changeRow3));
        changeRow4=TEMP2<croOptions.Obj.lb;
        TEMP2(find(changeRow4))=2*croOptions.Obj.lb(find(changeRow4))-TEMP2(find(changeRow4));
         
         fSwarm1= feval(f2eval,R7,Loninput);
         fSwarm2= feval(f2eval,R10,Loninput);
         fSwarm3= feval(f2eval,TEMP1,Loninput);
         fSwarm4= feval(f2eval,TEMP2,Loninput);
         A=[fSwarm1;fSwarm2;fSwarm3;fSwarm4];
         B=[R7;R10;TEMP1;TEMP2];
         [pe,x]=min(A);
         C=[En(f,1);En(h,1)];
         [temppemin,y]=min(C);
         if pe<temppemin 
             tag1=1;
         else
             tag1=0;
         end
          if(energytemp>pe)
              NumHit=NumHit+1;
              ke=energytemp-pe;
              if f<h
                Swarm(f,:)=B(x,:);
                En(f,1)=pe;
                En(f,2)=ke;
                if tag1==1
                    PBest(f,:)=B(x,:);
                end
                Swarm(h,:)=Swarm(popsize,:);
                PBest(h,:)=PBest(popsize,:);
                PBest(popsize,:)=[];
                En(h,1)=En(popsize,1);
                En(h,2)=En(popsize,2);
                Swarm(popsize,:)=[];
                En(popsize,:)=[];
              elseif f>h 
                Swarm(h,:)=B(x,:);
                En(h,1)=pe;
                En(h,2)=ke;
                if tag1==1
                    PBest(h,:)=B(x,:);
                end
                PBest(f,:)=PBest(popsize,:);
                PBest(popsize,:)=[];
                Swarm(f,:)=Swarm(popsize,:);
                En(f,1)=En(popsize,1);
                En(f,2)=En(popsize,2);
                Swarm(popsize,:)=[];
                En(popsize,:)=[];
              end
              popsize=popsize-1; 
              if  pe<MinPE1
                 if x==1
                     MinStruct=R7; 
                 end
                 if x==2
                     MinStruct=R10; 
                 end
                 if x==3
                     MinStruct=TEMP1; 
                 end
                 if x==4
                     MinStruct=TEMP2; 
                 end
                 MinPE1=pe; 
                 MinHit=NumHit; 
             end
          end
                  
 end
 function[MinPE1,En,Swarm]=interineffectivecollision(croOptions,f2eval,R2,R3,En,Swarm,MinPE1,Loninput)
       %global  fSwarm ; %x1 x2 x3 x4 x5 x6 x7 x8 x9 x10
       global MinHit;
       global MinStruct;
       global f;
       global h;
       global StepSize;
       global NumHit;
       global inter;
       global PBest;
       inter=inter+1;
     energytemp=En(f,1)+En(f,2)+En(h,1)+En(h,2);
     r=rand(1);
     if r>=croOptions.SParams.w_global
         tag=1;
     else
         tag=0;
     end
     for i=1:croOptions.Vars.Dim
        c=StepSize(:,i);
        C1(1,:)=normrnd(0,c);
        if(C1(1,:)<=croOptions.Obj.lb(:,i))
            C1(1,:)= croOptions.Obj.lb(:,i);
        elseif(C1(1,:)>=croOptions.Obj.ub(:,i))
            C1(1,:)=croOptions.Obj.ub(:,i);
        end
        C2(1,:)=normrnd(0,c);
        if(C2(1,:)<=croOptions.Obj.lb(:,i))
            C2(1,:)= croOptions.Obj.lb(:,i);
        elseif(C2(1,:)>=croOptions.Obj.ub(:,i))
            C2(1,:)=croOptions.Obj.ub(:,i);
        end
     end
        Q1 = rand(1, croOptions.Vars.Dim);
        W1 = rand(1, croOptions.Vars.Dim);
        R8 = R2+C1+tag*(croOptions.PSO.c1*Q1.*(PBest(f,:)-R2) + croOptions.PSO.c2*W1.*(MinStruct-R2));
        Q2 = rand(1, croOptions.Vars.Dim);
        W2 = rand(1, croOptions.Vars.Dim);
        R9 = R3+C2+tag*(croOptions.PSO.c1*Q2.*(PBest(h,:)-R3)+croOptions.PSO.c2*W2.*(MinStruct-R3));
        changeRow1=R8>croOptions.Obj.ub;
        R8(find(changeRow1))=2*croOptions.Obj.ub(find(changeRow1))-R8(find(changeRow1));
        changeRow2=R8<croOptions.Obj.lb;
        R8(find(changeRow2))=2*croOptions.Obj.lb(find(changeRow2))-R8(find(changeRow2));
        changeRow3=R9>croOptions.Obj.ub;
        R9(find(changeRow3))=2*croOptions.Obj.ub(find(changeRow3))-R9(find(changeRow3));
        changeRow4=R9<croOptions.Obj.lb;
        R9(find(changeRow4))=2*croOptions.Obj.lb(find(changeRow4))-R9(find(changeRow4));
        
    fSwarm5= feval(f2eval,R8,Loninput);
    fSwarm6= feval(f2eval,R9,Loninput);
    pe3=fSwarm5; 
    pe4=fSwarm6;
    if pe3<En(f,1)
        PBest(f,:)=R8;
    end
    if pe4<En(h,1)
      PBest(h,:)=R9;
    end
    n=pe3+pe4;
     if(energytemp-n>0)
         NumHit=NumHit+1;
         a=rand(1);
         ke3=(energytemp-n)*a;
         ke4=(energytemp-n)*(1-a);
         Swarm(f,:)=R8;
         En(f,1)=pe3;
         En(f,2)=ke3;
         Swarm(h,:)=R9;
         En(h,1)=pe4;
         En(h,2)=ke4;    
     end
     if((pe3<MinPE1)&&(pe3<=pe4))
         MinPE1=pe3;
         MinStruct=R8;
         MinHit=NumHit; 
      elseif(pe4<MinPE1)
         MinHit=NumHit; 
         MinPE1=pe4;
         MinStruct=R9;
      end

 end
 function[MinPE1,En,Swarm]=final(croOptions,f2eval,Swarm,En,MinPE1,Loninput)
      global MinHit;
    global MinStruct;
    global popsize;
    while(popsize>=2)
        tag=0;
        k=0;
        NewSwarm(1,:)=MinStruct;
        PE1(1,1)=MinPE1;
        KE1(1,1)=MinPE1;
        for i=1:popsize-1
                tempEn=En;
                tempSwarm=Swarm;
            for j=(i+1):popsize
                TEMP1=tempSwarm(i,:);
                TEMP2=tempSwarm(j,:);
                for n=1:croOptions.Vars.Dim
                    a=rand(1);
                    if(a>0.5)
                        R1(1,n)=TEMP1(1,n);
                        R2(1,n)=TEMP2(1,n);
                    else
                        R1(1,n)=TEMP2(1,n);
                        R2(1,n)=TEMP1(1,n);
                    end
                end
                fSwarm1= feval(f2eval,R1,Loninput);
                fSwarm2= feval(f2eval,R2,Loninput);
                A=[fSwarm1;fSwarm2];
                B=[TEMP1;TEMP2];
                [minpe,x]=min(A);
                if minpe<MinPE1
                    k=k+1;
                    MinPE1=minpe;
                    MinStruct=B(x,:);
                    NewSwarm(k,:)=MinStruct;
                    PE1(k,:)=minpe;
                    KE1(k,:)=tempEn(i,2);
                    MinHit=MinHit+1;
                    tag=1;
                end
            end
        end
        En=[PE1,KE1];
        Swarm=NewSwarm;
        for i=1:k
           NewSwarm(k-i+1,:)=[]; 
        end
        [row,column]=size(Swarm);
        popsize=row;
    end
 end


