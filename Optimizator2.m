classdef Optimizator2
    
   properties
      optimization_factor;
      optimization_factor_last;
      num_loops;
      range_optimization_search_sigma; %percent value of sigma for search of optimum value
       num_of_variable_to_optimize; %modulus and angle of SIGMA
      last_evaluation_index;
      evaluation_index;
      internal_step;
      last_perturbed_variable;
      grad;
      alpha_min;
      alpha_max;
      perturbation;
     
   end
   
   methods
       function obj = Optimizator2(numLoops)
          obj.num_loops = numLoops;
          
         
          obj.optimization_factor(1:2*obj.num_loops) = 100*ones(2*obj.num_loops,1); %percent of modulus and angle
           obj.optimization_factor_last(1:2*obj.num_loops) = 100*ones(2*obj.num_loops,1); %percent of modulus and angle
          
        obj.range_optimization_search_sigma = 90; %percent value of sigma for search of optimum value
        obj.perturbation = 0.02;
       obj.num_of_variable_to_optimize = 2*obj.num_loops; %modulus and angle of SIGMA
       obj.last_evaluation_index = 0;
       obj.internal_step = 1;
       obj.last_perturbed_variable=0;
       obj.grad(1:2*obj.num_loops) = 100*zeros(2*obj.num_loops,1);
       obj.alpha_min = 1000;
       obj.alpha_max = 0;
       end
       
       function obj = step(obj, eval_index)
            
            obj.evaluation_index = eval_index;
            
             if obj.internal_step==1
           obj.last_evaluation_index = obj.evaluation_index;
           end
            
            obj.optimization_factor = obj.optimization_factor_last;
           
            if obj.internal_step>0 && obj.internal_step <= obj.num_of_variable_to_optimize
                obj.optimization_factor(obj.internal_step) = obj.optimization_factor(obj.internal_step) * (1 + obj.perturbation);
            end
            
            if obj.last_perturbed_variable > 0 && obj.last_perturbed_variable <= obj.num_of_variable_to_optimize
               obj.grad(obj.last_perturbed_variable) = (obj.evaluation_index - obj.last_evaluation_index) / (obj.perturbation * obj.optimization_factor(obj.last_perturbed_variable));
               % obj.grad(obj.last_perturbed_variable) = (obj.evaluation_index - obj.last_evaluation_index);
               obj.alpha_min = min([obj.alpha_min abs(obj.evaluation_index - obj.last_evaluation_index)]);
               obj.alpha_max = max([obj.alpha_max abs(obj.evaluation_index - obj.last_evaluation_index)]);
            end
            
            
            obj.last_perturbed_variable = obj.internal_step;
                    
           
          
            if obj.internal_step > obj.num_of_variable_to_optimize
               if (obj.alpha_min==0) 
                   obj.alpha_min=0.1;
               end
               if (obj.alpha_max==0) 
                   obj.alpha_max=0.1;
               end
               %alpha = (abs(0.01/obj.alpha_min) + abs(0.08/obj.alpha_max))/5;
               alpha = 5*(abs(obj.alpha_min) + abs(obj.alpha_max));
               
               for u=1 : obj.num_of_variable_to_optimize
                   var = abs( alpha * obj.grad(u));
                   if var <0.005
                       var = 0.005;
                   end
                   if var >0.016
                       var = 0.16;
                   end
                   %var = 0.02;
                   var = var * sign(alpha * obj.grad(u)) * -1;
                       
                   obj.optimization_factor(u) = obj.optimization_factor(u) + obj.optimization_factor(u) *var;
                         if obj.optimization_factor(u) > 100+obj.range_optimization_search_sigma
                                obj.optimization_factor(u) = 100+obj.range_optimization_search_sigma;
                         elseif obj.optimization_factor(u) < 100-obj.range_optimization_search_sigma 
                                obj.optimization_factor(u) = 100-obj.range_optimization_search_sigma;
                         end

                         if obj.optimization_factor(u) < 0 
                                obj.optimization_factor(u) = 0;
                         end
                   
               end
               obj.optimization_factor_last = obj.optimization_factor;
               obj.alpha_min = 1000;
               obj.alpha_max = 0;
           end
                 
           obj.internal_step = obj.internal_step + 1;
           
           if obj.internal_step > obj.num_of_variable_to_optimize + 1
               obj.internal_step = 1;
           end   
          
           
                 
                 
       end
      
   end
   
end