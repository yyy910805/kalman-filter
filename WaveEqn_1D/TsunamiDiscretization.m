classdef TsunamiDiscretization < noname.Discretization
    properties
        name         = 'Shallow water with shore at left boundary and open water at right boundary'
        description  = 'A*u_t + B*u_x = G'
        order        %Order of accuracy
        
        D            %Discretization matrix including BC
        S            %Function handle, v_t = D*v + S(t)
        F            %Function handle, v_t = F(v,t) = D*v + S(t)
        v0           %Initial data
        diffOp       %Various relevant operators      
    end

    methods
        % p on primal grid. Here, p = q = flux.
        % v on dual grid. Here, v = h = wave height.
        % G: forcing function, cell array {Gp(x,t), Gv(x,t)} (pressure and velocity forcings)
        function obj = TsunamiDiscretization(m, order, lim, G, A, B, boundary_data_l, boundary_data_r)
            default_arg('boundary_data_r',[]);
            default_arg('boundary_data_l',[]);
            default_arg('B',{@(x)0*x, @(x)0*x+1; @(x)0*x+1, @(x)0*x});
            default_arg('A',{@(x)0*x+1, @(x)0*x; @(x)0*x, @(x)0*x+1});
            default_arg('G',[]);
            default_arg('lim',{0,1});
            default_arg('order',4);
            
            % Create grid and differential operators (diffOp) on grid
            [g_primal, g_dual] = grid.primalDual1D(m, lim);
            g = grid.Staggered1d(g_primal, g_dual);
            diffOp = scheme.Staggered1DAcousticsVariable(g, order, A, B);
            D = diffOp.D;

            % Forcing function G
            x_primal = diffOp.grid_primal.points();
            x_dual = diffOp.grid_dual.points();
            if ~isempty(G)
                Gt = @(t) [(1./A{1,1}(x_primal)).*G{1}(x_primal,t);...
                           (1./A{2,2}(x_dual)).*G{2}(x_dual,t)];
            else
                Gt = [];
            end

            % Set boundary conditions
            [closure_l, penalty_l] = diffOp.boundary_condition('l','p');
            [closure_r, penalty_r] = diffOp.boundary_condition('r','characteristic');
            D = D + closure_l + closure_r;
            obj.D = D;

            % Create data function S(t)
            if ~isempty(boundary_data_l)
                boundary_data_l = @(t) penalty_l*boundary_data_l(t);
            end
            if ~isempty(boundary_data_r)
                boundary_data_r = @(t) penalty_r*boundary_data_r(t);
            end
            data_funcs = {Gt, boundary_data_l, boundary_data_r};
            
            S = [];
            for i = 1:numel(data_funcs)
                data = data_funcs{i};
                if ~isempty(data)
                    if(isempty(S))
                        S = data;
                    else
                        S = @(t) S(t) + data(t);
                    end
                end
            end
            obj.S = S;        

            % Initial data
            x0 = 0.5;
            sigma = 0.05;
            p0_fun = @(x)0*x;
            v0_fun = @(x)exp(-((x-x0).^2)/sigma^2);
            obj.v0 = [p0_fun(x_primal); v0_fun(x_dual)];
            obj.diffOp = diffOp;
            
            if isempty(obj.S)
                obj.F = @(v,t) obj.D*v;
            else
                obj.F = @(v,t) obj.D*v + obj.S(t);
            end

        end
        % Prints some info about the discretisation
        function printInfo(obj)
            fprintf('Name: %s\n',obj.name);
            fprintf('Size: %d\n',obj.size());
        end

        % Return the number of DOF
        function n = size(obj)
            n = length(obj.v0);
        end

        % Returns a timestepper for integrating the discretisation in time
        %     method is a string that states which timestepping method should be used.
        %          The implementation should switch on the string and deliver
        %          the appropriate timestepper. It should also provide a default value.
        %     time_align is a time that the timesteps should align with so that for some
        %                integer number of timesteps we end up exactly on time_align
        function [ts, N] = getTimestepper(obj,method,time_align) 
            default_arg('method','rk4');
            default_arg('time_align',[]);
            switch method
                case 'rk4'
                    k = obj.getTimestep(method);

                    if ~isempty(time_align)
                        [k, N] = alignedTimestep(k, time_align);
                    end

                    t = 0;
                    ts = time.Rungekutta4proper(obj.F,k,t,obj.v0);
                otherwise
                    error('Timestepping method ''%s'' not supported',method);
            end
        end

        function k = getTimestep(obj, method, cfl)
            k = 0.2*time.rk4.get_rk4_time_step(max(abs(eig(full(obj.D)))));
        end

        function r = getTimeSnapshot(obj, ts)
            if ts == 0
                r.t = 0;
                r.v = obj.v0;
                return
            end
            r.t = ts.t;
            r.v = ts.getV();
        end

        % Sets up movie recording to a given file.
        %     saveFrame is a function_handle with no inputs that records the current state
        %               as a frame in the moive.
        function saveFrame = setupMov(obj, file)
            error('not implemented');
        end

        % Sets up a plot of the discretisation
        %     update is a function_handle accepting a timestepper that updates the plot to the
        %            state of the timestepper
        function [update,figure_handle] = setupPlot(obj, type)
            x_primal = obj.diffOp.grid_primal.points();
            x_dual = obj.diffOp.grid_dual.points();

            e_primal = obj.diffOp.e_primal;
            e_dual = obj.diffOp.e_dual;

            figure_handle = figure();

            line_p = plot(x_primal, e_primal'*obj.v0, 'r-o');
            hold on
            line_v = plot(x_dual, e_dual'*obj.v0, 'b-o');
            ylim([-1, 1])
            legend('q','h')
            figure(figure_handle)

            a = gca;

            function update_fun(r, e_primal, e_dual)
                t = r.t;
                v = r.v;
                if ishandle(a)
                    title(a,sprintf('T = %.3f',t))
                    line_p.YData = e_primal'*v;
                    line_v.YData = e_dual'*v;
                end
            end
            update = @(r)update_fun(r, e_primal, e_dual);
        end

        % Compare the gridfunctions to the analytical functions gp(x) and gv(x) in the sbp-norm (weighted l2-norm).
        function [e, e_vec] = compareSolutionsAnalytical(obj, u, gp, gv)

            x_primal = obj.diffOp.grid_primal.points();
            x_dual = obj.diffOp.grid_dual.points();

            p_exact = gp(x_primal);
            v_exact = gv(x_dual);
            u_exact = [p_exact; v_exact];

            H = obj.diffOp.H;
            e_vec = u-u_exact;
            e = sqrt(e_vec'*H*e_vec);
        end

    end

    methods(Static)
        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(u, v)
            error('not implemented');
        end
    end
end