%%% THIS CODE IS THE MAIN OPTIMIZATION LOOP FOR OBJECTIVE AND CONSTRANTS
%%% Date: 06042020
function [doptim_curr,foptim,eflag]=solve_inneroptimzation(objfunc, nonlcon,NDV,lb_design,ub_design,.....
    A,b,Aeq,beq,solver_choice,opts1,opts2_xinit)

if strcmp(solver_choice,'ga')
    %% carry out optimization using the surrogates to find the optimum!!!
    %% use GA for optimization
    options = optimoptions('ga','UseParallel',true,'MaxGenerations',opts1,'ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf);
    [doptim_curr,foptim,eflag]=ga(objfunc,NDV,A,b,Aeq,beq,lb_design,ub_design,nonlcon,options);
elseif strcmp(solver_choice,'ms')
    %% multistart global optimization!!!!
    if isempty(opts2_xinit)
        xm_initial=lb_design+(ub_design-lb_design).*rand(1,NDV);
    else
        xm_initial=opts2_xinit ;    %% can update this iteratively based on the previous optimum!!!!!!!
    end
    npts_ms=opts1;  %number of points for multistart
    
%     [objfunc(xm_initial) nonlcon(xm_initial)]
        opts = optimoptions(@fmincon,'Algorithm','sqp','ConstraintTolerance',1e-6,'UseParallel',true,'display','off',...
            'MaxIterations',2.5e4,'MaxFunctionEvaluations',5e5 );
%     opts = optimoptions(@fmincon,'Algorithm','interior-point','ConstraintTolerance',1e-6,'UseParallel',true);
    problem = createOptimProblem('fmincon','x0',xm_initial,'objective',objfunc,'lb',lb_design,'ub',ub_design,....
        'Aineq',A,'bineq',b,'Aeq',Aeq,'beq',beq,'lb',lb_design,'ub',ub_design,'nonlcon',nonlcon,'options',opts);
    %     ms = MultiStart('FunctionTolerance',2e-4,'XTolerance',5e-3,...
    %         'StartPointsToRun','bounds-ineqs')
%     ms = MultiStart('FunctionTolerance',1e-6,'XTolerance',1e-6,'StartPointsToRun','bounds');
    ms = MultiStart('FunctionTolerance',1e-6,'XTolerance',1e-6);
    [doptim_curr,foptim,eflag] = run(ms,problem,npts_ms);
    
end