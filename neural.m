%net = feedforwardnet(2,'trainlm');



%teacher = rand(100,100);
%output = rand(1,100);
%net = train(net,teacher, output);
%


% learning_rate = 0.2;
% dataset = randn(1,1000);
% dataset1 = randn(1,1000);
% teacher_w1 = randn(1000,1);
% student_w1 = randn(1000,1);
% teacher_w2 = randn(1000,1);
% student_w2 = randn(1000,1);


%your_result = [];
% w_dif1 = sum(teacher_w1)- sum(student_w1)
%%
student_weights = [];
student_scalar = [];
teacher_scalar = [];
R_in = [];
Q_ik = [];
mu = 200000;
si = 0;
tou = 0;
samples= 100;
K = 2;
teacher_weights = randn(K,samples);
teacher_weights_ortho = GramSchmidt(teacher_weights');
teacher_weights = teacher_weights_ortho';

student_weights = randn(K,samples);

%for jj = 1:K
 %  teacher_weights(K,:) =  (teacher_weights(K,:))./(norm(teacher_weights(K,:)));
 %  student_weights(K,:) =  (student_weights(K,:))./(norm(student_weights(K,:)));
    
%end

for i = 1:mu

    dataset = randn(samples,1);
    for j = 1:K
   %i
    
   % teacher_1 = dot(teacher_w1, dataset);
    %student_1 = dot(student_w1, dataset);
    teacher_scalar = dot(teacher_weights(j,:), dataset);
    student_scalar(j) = dot(student_weights(j,:), dataset) ;

    activation_teacher  = erf(teacher_scalar/sqrt(2));
    activation_student  = erf(student_scalar(j)/sqrt(2));
    %activation_final_teacher = activation_final_teacher + activation_teacher;
    si = si + activation_student;
    tou = tou + activation_teacher;
    
    end
    
    for j = 1:K
        if mod(i,100) == 0
            R_in(i/100,j) = dot(student_weights(1,:), teacher_weights(j,:));
            Q_ik(i/100,j) = dot(student_weights(1,:), student_weights(j,:));
        end
        student_s = student_scalar(j);
        
        gradeint_epsilon_studentweights = (si - tou) * sqrt(2/pi)*exp(-(student_s*student_s)/2)* dataset;
        student_weights(j,:) =  student_weights(j,:) - ((1.5/samples)* gradeint_epsilon_studentweights');
         
        
    end
    
    si = 0;
    tou = 0;
    
end

x = 0:mu/100-1;

plot(x,R_in(:,1));
hold on
plot(x,R_in(:,2));
plot(x,Q_ik(:,1));
plot(x,Q_ik(:,2));
legend('R','R','Q','Q');
xlabel('Time steps') 
ylabel('Order parameter') 

%%
for j = 1:2
    %j
    dataset = randn(1,100); % this is definetly random all the time.
    teacher_weight = randn(100,1);
    student_weight = randn(100,1);
    
    for i = 1:mu

   %i
    
   % teacher_1 = dot(teacher_w1, dataset);
    %student_1 = dot(student_w1, dataset);
    teacher_scalar = dot(teacher_weight, dataset);
    student_scalar = dot(student_weight, dataset) ;

    activation_teacher  = erf(teacher_scalar/sqrt(2));
    activation_student  = erf(student_scalar/sqrt(2));
    %activation_final_teacher = activation_final_teacher + activation_teacher;
    si = si + activation_student;
    tou = tou + activation_teacher;
    epsilon = (1/2)*(si-tou)^2;
    gradeint_epsilon_studentweights = (si - tou) * sqrt(2/pi)*exp(-(student_scalar*student_scalar)/2).* dataset;
    gradeint_epsilon_studentweights;
    updated_student_weight = student_weight + ((0.5/100).* gradeint_epsilon_studentweights);
    %student_weight_inloop = student_weight; % initialising the student weight in the loop for 1st time to starting value of randn function.
    %student_weight_update = student_weight 
    %activation_final_teacher = activation_final_teacher + activation_student; 
    % 
%    g_teacher_2 = erf(teacher_2/sqrt(2));
%    g_student_2 = erf(student_2/sqrt(2));


    
   % gradient = sqrt(2/pi)*exp(-(student_1*student_1)/2);
   % lets try and calculate gradeint
   %%%%%%%%%%
   % TODO
   % how to store the value of si and tou in a matrix or a way that the
   % final equation of gradient becomes understandable.
   %%%%%%%%%%%
   
    %gen_error = gradient*0.5*(g_teacher_1 + g_teacher_2 - g_student_1 - g_student_2)^2;
    % what if I calculate the generilzation error outside the loop.
   % student_w1 = student_w1 - gen_error*learning_rate;
    %student_1 = dot(student_w1, dataset)/1000 ;
    %err_student_1 = erf(student_1/sqrt(2));
    
   % gradient = sqrt(2/pi)*exp(-(student_2*student_2)/2)
    
    %gen_error = gradient*0.5*(g_teacher_1 + g_teacher_2 - g_student_1 - g_student_2)^2;
    
    %student_w2 = student_w2 - gen_error*learning_rate;
    
    %y = linspace(1,100,100);
   % scatter(y, gradeint_epsilon_studentweights)
    %plot([y(end) y], [gradeint_epsilon_studentweights(end) gradeint_epsilon_studentweights], 'r-');
    %hold on
    end
    % try to append the pdated weights 
    %updated_trial_student_weights = [updated_trial_student_weights;updated_student_weight];
end
%Something
updated_trial_student_weights;
% extracting the last column of updated_trial_student_weights and then
% plotting against teacher weight.
y = linspace(1,100,100);
%scatter(y, updated_trial_student_weights[:,end]);
 
%w_dif2 = sum(teacher_w1)- sum(student_w1);
