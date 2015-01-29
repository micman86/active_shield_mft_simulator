
close all
clear variables


%% PARAMETERS

%substation_mvlv_4loops_parameters;
example_1loop_parameters;

num_points = size(alpha_physical,1);
num_loops = size(beta_physical,2);

optimizator = Optimizator2(num_loops);






%% SIMULATION

%% initial values


num_of_variable_to_optimize = 2*num_loops; %modulus and angle of SIGMA

source_current_index(1) = 1;
sensed_current_index(1) = 1;
shielding_current(1,1:num_loops) = zeros(num_loops,1);
set_shielding_current(1,1:num_loops) = zeros(num_loops,1);

source_current(1) = source_current_array(source_current_index,2);
sensed_current(1) = source_current_array(source_current_index,2);
index = 1;
time_array(1) = 0;


sigma = sigma_initial; 
optimization_factor(1,1:2*num_loops) = 100*ones(2*num_loops,1); %percent of modulus and angle 
evaluation_index(1) = 0;


%for opt_index = 100 - range_optimization_search_sigma: 1: 100 + range_optimization_search_sigma


% %sigma_array = repmat(sigma_class(), 100 + range_optimization_search_sigma);
% for var=1 : num_of_variable_to_optimize
%     for opt_index=100 - range_optimization_search_sigma: 100 + range_optimization_search_sigma
%         sigma_array(var,opt_index) = sigma_class(var,opt_index);
%     end
% end



%% simulation
for time = shielding_time_step: shielding_time_step: simulation_duration
    
    index = index + 1;
    time_array(index) = time;
    
    %% Iso current SET
    if size(source_current_array,1) == source_current_index
        source_current(index) = source_current(index-1);
    else
        if source_current_array(source_current_index+1, 1) < time
            source_current_index = source_current_index + 1;
            source_current(index) = source_current_array(source_current_index,2);
        else 
            source_current(index) = source_current(index-1);
        end
    end
    
    
     if size(source_current_array,1) == sensed_current_index
        sensed_current(index) = sensed_current(index-1);
    else
        if source_current_array(sensed_current_index+1, 1) < time - shielding_time_step
            sensed_current_index = sensed_current_index + 1;
            sensed_current(index) = source_current_array(sensed_current_index,2);
        else 
            sensed_current(index) = sensed_current(index-1);
        end
    end
    
   
    
    %% alghoritm simulation
    
    shielding_current(index,:) = Ish_driver_factor*set_shielding_current(index-1,:);
    
    %detected_B(index) = B_probe_factor * (alpha_physical*sensed_current(index) + beta_physical*shielding_current(index));
    
    for point=1:num_points
    detected_B(point,index) = B_probe_factor * (alpha_physical(point)*sensed_current(index) + sum(beta_physical(point,:).*shielding_current(index,:)));
    real_B(point,index) = (alpha_physical(point)*source_current(index) + sum(beta_physical(point,:).*shielding_current(index,:)));
    
    
    for i_axis=1:3
        B_rms(point,index,i_axis) = alpha_physical_3d(point,i_axis)*source_current(index);
        for loop=1:num_loops
            B_rms(point,index,i_axis) = B_rms(point,index,i_axis) + beta_physical_3d(point,loop,i_axis)*shielding_current(index,loop);
        end
     end
    
    B_tot(point,index) = abs(sqrt(sum(abs(B_rms(point,index,1:3)).^2)));
    
    end
    detected_Iso(index) = (detected_B(1,index) - I_probe_factor*sum(shielding_current(index).*beta(1,:)))/alpha(1);
    
   
    
    
    %% optimization
    
      if enable_optimization == 1
          
          evaluation_index(index) = sum(weights_eval .* B_tot(:,index));
                
            if (index > 50 && isnan(detected_B(index)) == false)
                
                optimizator = optimizator.step(evaluation_index(index));
                optimization_factor(index,:) = optimizator.optimization_factor;
                
            else
                optimization_factor(index,:) = optimization_factor(index-1,:);
            end 


                    %% set Ish
            for loop=1:num_loops
                     sigma(loop) = (abs(sigma_initial(loop)) * optimization_factor(index,1+2*(loop-1)) / 100)*exp(1i*angle(sigma_initial(loop)) * optimization_factor(index,2+2*(loop-1)) / 100);
            end  
        end


         set_shielding_current(index,:) = sigma .* detected_Iso(index);

    
    
end



    
 
 %% matlab plots and outputs

% %h = plot(time_array,abs(source_current),time_array,abs(detected_Iso),time_array,abs(real_B),time_array,abs(detected_B),time_array,abs(set_shielding_current),time_array,evaluation_index,time_array,optimization_factor(:,1),time_array,optimization_factor(:,2),'linewidth',3);
% h = plot(time_array,abs(evaluation_index),time_array,abs(detected_B),time_array,abs(set_shielding_current),time_array,optimization_factor(:,1),time_array,optimization_factor(:,2),'linewidth',3);
% legend(h,'actual B (uT)','detected B (uT)', 'set Ish (A)','Perturbation modulus (%)','Perturbation phase (%)')
% xlabel('Time (s)','Fontsize',18);
% set(gca,'FontSize',18)
% grid on
% set(gcf,'color',[1 1 1]);

% h2 = plot(time_array,detected_Iso,time_array,evaluation_index,time_array,B_tot(1,:),time_array,B_tot(2,:),time_array,B_tot(3,:),time_array,optimization_factor(:,1),time_array,optimization_factor(:,2),time_array,optimization_factor(:,5),time_array,optimization_factor(:,6),'linewidth',3);
% legend(h2,'Iso','eval idx','B 1 (uT)','B 2 (uT)', 'B 3 (uT)','Perturbation modulus 1 (%)','Perturbation phase 1 (%)','Perturbation modulus 3 (%)','Perturbation phase 3 (%)')
% xlabel('Time (s)','Fontsize',18);
% set(gca,'FontSize',18)
% grid on
% set(gcf,'color',[1 1 1]);




% h2 = plot(time_array,evaluation_index,time_array,B_tot(1,:),time_array,B_tot(2,:),time_array,B_tot(3,:),time_array,B_tot(4,:),time_array,optimization_factor(:,1),time_array,optimization_factor(:,2),time_array,optimization_factor(:,5),time_array,optimization_factor(:,6),'linewidth',3);
% legend(h2,'eval idx','B 1 (uT)','B 2 (uT)', 'B 3 (uT)', 'B 4 (uT)','Perturbation modulus 1 (%)','Perturbation phase 1 (%)','Perturbation modulus 3 (%)','Perturbation phase 3 (%)')
% xlabel('Time (s)','Fontsize',18);
% set(gca,'FontSize',18)
% grid on
% set(gcf,'color',[1 1 1]);



h2 = plot(time_array,evaluation_index,time_array,B_tot(1,:),time_array,optimization_factor(:,1),time_array,optimization_factor(:,2),'linewidth',3);
legend(h2,'eval idx','B 1 (uT)','Perturbation modulus 1 (%)','Perturbation phase 1 (%)')
xlabel('Time (s)','Fontsize',18);
set(gca,'FontSize',18)
grid on
set(gcf,'color',[1 1 1]);


% h3 = figure; set(gcf,'color',[1 1 1])
% plot(time_array,abs(source_current),time_array,abs(set_shielding_current),time_array,abs(real_B),'linewidth',3);
% legend2 = legend('actual Iso (A)', 'set Ish (A)', 'actual B (\muT)');
% xlabel('Time (s)','Fontsize',18);

% title('Active magnetic shielding with online optimization','Fontsize',18);
set(gca,'FontSize',18)
grid on

sprintf('shielding currents')
shielding_current(index,:)

%sigma
sprintf('B tot')
B_tot(:,index)

