%% The script plots the restarts method.
% LR can be changed to observe different order parameters respectively
% initialization of variables and vectors
student_scalar = [];
teacher_scalar = [];
gen_error = [];
R_in = [];
Q_ik = [];
mu = 200000;
si = 0;
tou = 0;
samples= 1000; % No of examples or samples or Dimensions
% hyperparameter gamma
gamma = 0.1
learning_rate = 1.5;
beta = 0.9
step = 10;
K=2; % no of hidden units
epsilon = 0;
cycles = 0;
rng(5782)
teacher_weights = randn(K,samples)./sqrt(samples); % teacher weight vector
teacher_weights_ortho = GramSchmidt(teacher_weights'); %conversion toIsotropic teachers
teacher_weights = teacher_weights_ortho';
student_weights = randn(K,samples)./sqrt(samples);

velocities = 0;
for i = 1:mu
    
    % same array of input for the different weight vectors.
    % input to be generated for all the runs?
    
    
    %Implementation of SGDR using cosine annealing
    if  learning_rate > 0.1
      %  learning_rate = learning_rate - 1/(samples*10)
       learning_rate = cos(epsilon)*1.5;
       epsilon = epsilon + step/samples;
    else
        learning_rate = 1.5;
        epsilon = 0;
        step = step/2;
    end
    
    % creating datavector for every run using a random seed.
    rng(rngSeed)
    dataset = randn(samples,1);
    rngSeed = rngSeed + 1;
    for j = 1:K
    % dot product of teacher and student weights with the input dataset.
    teacher_scalar = dot(teacher_weights(j,:), dataset);
    student_scalar(j) = dot(student_weights(j,:), dataset) ;
    
    % sigmoid activation function
    activation_teacher  = erf(teacher_scalar/sqrt(2));
    activation_student  = erf(student_scalar(j)/sqrt(2));
    si = si + activation_student;
    tou = tou + activation_teacher;
    
    end
    % both layers getting updated simultaneously for mu times.
    for j = 1:K
        if mod(i,samples) == 0
            % R = 1st student dot product with both teacher weight vetor
            R_in(i/samples,j) = dot(student_weights(1,:), teacher_weights(j,:));
            % Q = 1st student dot product with itself and other student
            % weight vetor
            Q_ik(i/samples,j) = dot(student_weights(1,:), student_weights(j,:));
            % gen_error = vector for number of steps. generalisation error
            % is for both the weight vectors together in one vector.
            gen_error(i/samples) = generalizationerror(student_weights,teacher_weights,K);
        end
        student_s = student_scalar(j);
        
        % Calculating the gradient 
        gradient_epsilon_studentweights = (si - tou) * sqrt(2/pi)*exp(-(student_s*student_s)/2)* dataset;
         velocities = beta * velocities + (1-beta)* gradient_epsilon_studentweights';
         student_weights(j,:) =  student_weights(j,:) - learning_rate/samples*velocities;
        
    end
    
    si = 0;
    tou = 0;
    
end


%plotting Order Parameters
x = 0:mu/samples-1;
figure(1)
plot(x,R_in(:,1));
hold on

plot(x,R_in(:,2));
plot(x,Q_ik(:,1));
plot(x,Q_ik(:,2));
xlabel('alpha')
ylabel('Order-Parameters')
legend('R','R','Q','Q')
gen_error;

%Plotting generalization error

figure(2)
plot(x, gen_error)
xlabel('alpha')
ylabel('Generalization error')
legend('Restarts')
hold on


