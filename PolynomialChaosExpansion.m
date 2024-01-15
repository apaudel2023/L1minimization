classdef PolynomialChaosExpansion
    properties
        func
        dim
        max_order
        a
        b
        result_file_name
    end

    methods

        function obj = PolynomialChaosExpansion(func, dim, max_order, a, b, result_file_name)
            obj.func = func;
            obj.dim = dim;
            obj.max_order = max_order;
            obj.a = a;
            obj.b = b;
            obj.result_file_name = result_file_name;
        end

        % Generate Samples from the distribution
        function [standardized, original] = generateSamples(obj, nsamp)
            basis = lhsdesign(nsamp, obj.dim);

            % Transform LHS basis to range [-1, 1] for standardized samples
            standardized = 2 * basis - 1; % This transforms [0, 1] to [-1, 1]

            % Transform LHS basis to specified bounds for original samples
            original = zeros(size(basis));
            for i = 1:obj.dim
                % Check if bounds are valid
                if obj.a(i) >= obj.b(i)
                    error('Invalid bounds: Lower bound must be less than upper bound.');
                end

                % Transform basis to the range specified by a and b
                original(:, i) = basis(:, i) * (obj.b(i) - obj.a(i)) + obj.a(i);
            end
        end

        % Build Amatrix for solving pce coefficient
        % uses multiindices and basis samples
        % basis samples example: U[-1,1] and N[0,1]
        function amatrix = computeAmatrix(obj, multiIndices, basisSamp)
            DIM = size(multiIndices,2);
            amatrix = zeros(size(basisSamp,1),size(multiIndices,1));
            for rowitk=1:size(basisSamp,1)
                for colitk=1:size(multiIndices,1)
                    ii = 1;
                    amatrix(rowitk,colitk) = 1;
                    while (ii<=DIM)
                        amatrix(rowitk,colitk)=amatrix(rowitk,colitk)*LEGENPOLY(multiIndices(colitk,ii),basisSamp(rowitk,ii));
                        ii = ii+1;
                    end
                end
            end
        end

        % Get MCS reponse 
        % samples must be in original domain
        function response_mcs = mcsResponse(obj, x_mcs_original)
            response_mcs = obj.func(x_mcs_original);
        end
 
        % ==========================================
        % Define different PCE solutions schemes
        % L2 based, L1 based, L1/L2 combination

        % PCE with L1-minimization using optimization: 
        function [pce_coeff_final, best_lambda] = solveL1Minimization(obj, Amatrix, Amatrix_test, response_build_pce, mcs_response)
            % Discretize the regularization parameter
            lamdamax = 1;
            lamdamin = 1e-6;
            nlamda=5;  %total number of increments for the lamda parameter!!!!
            lambdaArr = sort(logspace(log10(lamdamin),log10(lamdamax),nlamda),'descend')';

            % set up the optimization problem
            lb_design = [];
            ub_design = [];
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            nonlcon = [];
            solver_choice = 'ms';
            NDV = size(Amatrix, 2); % no of design variable

            maxGenGA = 5e2;
            msNum = 10;

            if solver_choice == 'ga'
                opts1 = maxGenGA;
                opts2_xinit = [];
            elseif solver_choice == 'ms'
                opts1 = msNum;
                a0 = zeros(1,size(Amatrix, 2));
                opts2_xinit = a0';
            end

            % get pce coeffcients for each lambda
            for iit = 1:length(lambdaArr)
                lambda = lambdaArr(iit);
           
                % minimization objective
                objective_func=@(inputs) norm(response_build_pce-Amatrix*(inputs)) / 2 + lambda*sum(abs(inputs));

                [optX, foptim, exitflag] = solve_inneroptimzation(objective_func, nonlcon, NDV,lb_design, ub_design,...
                    A, b, Aeq, beq, solver_choice, opts1, opts2_xinit);
                pce_coeff_l1(:, iit) = optX;
                pce_response_test_l1(:,iit) = Amatrix_test * pce_coeff_l1(:,iit);
                kld_l1(iit,1) = calculateKLDivergence(mcs_response, pce_response_test_l1(:,iit));
            end

            % select the best solution based on kld error
            [sorted_array, sort_order] = sort(kld_l1);
            best_lambda = lambdaArr(sort_order(1));
            pce_coeff_final = pce_coeff_l1(:, sort_order(1));
            
        end
        
        function solveL1Analytical(obj)
            % Discretize the regularization parameter
            lamdamax = 1;
            lamdamin = 1e-6;
            nlamda=5;  %total number of increments for the lamda parameter!!!!
            lambdaArr = sort(logspace(log10(lamdamin),log10(lamdamax),nlamda),'descend')';


        end

        % PCE with L2 minimization (Least-Squares)
        function x = solveLeastSquares(obj, Amat, y)
            % Solve the least-squares
            x = linsolve(Amat, y);
        end


        % Get PCE response based on the obtained coefficients 
        function pce_response_test = getPCEresponse(obj, Amatrix, coeff)
            pce_response_test = Amatrix * coeff;
        end


        % ============== Main Loop: Run Everything ============
        function runPCE(obj, nsamp_pce, nsamp_mcs)
            fprintf('PCE Initiated ...\n')

            % for building pce
            [basis_samp_build_pce, org_samp_build_pce] = obj.generateSamples(nsamp_pce);
            response_build_pce = obj.func(org_samp_build_pce);
            Amatrix = [];

            % for testing pce
            Amatrix_test = [];
            [basis_samp_test_pce, org_samp_test_pce] = obj.generateSamples(nsamp_mcs);
            mcs_response = obj.mcsResponse(org_samp_test_pce);

            for order = 1:obj.max_order
                fprintf('... order - %d\n', order)

                multiind_current = generateMultiIndices(obj.dim, order);

                Amatrix_current = obj.computeAmatrix(multiind_current, basis_samp_build_pce);
                Amatrix = [Amatrix Amatrix_current];

                Amatrix_current_test = obj.computeAmatrix(multiind_current, basis_samp_test_pce);
                Amatrix_test = [Amatrix_test Amatrix_current_test];

                % ====================
                % here you can change the type of solver to use for pce solutions
                pce_coeff = obj.solveLeastSquares(Amatrix, response_build_pce);

                % [pce_coeff, lambda]= obj.solveL1Minimization(Amatrix, Amatrix_test, response_build_pce, mcs_response);

                % ====================

                pce_response = obj.getPCEresponse(Amatrix_test, pce_coeff);
                kld = calculateKLDivergence(mcs_response, pce_response);

                PCE.response_pce(:, order) = pce_response;
                PCE.coeff{order, 1} = pce_coeff;
                PCE.kld(order, 1) = kld;
                if exist('lambda', 'var')
                    PCE.best_lambda = lambda;
                end

            end
            PCE.Amatrix_test = Amatrix_test;
            PCE.response_mcs = mcs_response;
            PCE.samples_basis = basis_samp_test_pce;
            PCE.samples_original = org_samp_test_pce;
            save(obj.result_file_name,'PCE')

            fprintf('Maximum order achieved ... Done !!! \n')
        end

    end
end