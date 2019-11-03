%% The script plots the orderparameter and generalization error for Momentum method 
% variable and vector definitions
student_scalar = [];
teacher_scalar = [];
gen_error = [];
R_in = [];
Q_ik = [];
mu = 140000;
si = 0;
tou = 0;
samples= 1000;% No of examples or samples or Dimensions
learning_rate = 1.5; % LR can be altered to observe differet OP

beta = 0.99 

K=2 % 2 hidden units used,
rngSeed = 999
teacher_weights = randn(K,samples)./sqrt(samples);% teacher weight vector
teacher_weights_ortho = GramSchmidt(teacher_weights');%conversion to Isotropic teachers
teacher_weights = teacher_weights_ortho';
student_weights = randn(K,samples)./sqrt(samples);

velocities = 0;
epsilon = 2*pi

for i = 1:mu

    
    if  beta > 0
      beta = cos(epsilon)*0.9;
      epsilon = epsilon - pi/10000;
    else
        beta = 0;
    end
    % creating input vector(si) for every run using a random seed.
    rng(rngSeed);
    dataset = randn(samples,1);
    rngSeed = rngSeed + 1;
    for j = 1:K
     % dot product of teacher and student weights with the input dataset.
    teacher_scalar = dot(teacher_weights(j,:), dataset);
    student_scalar(j) = dot(student_weights(j,:), dataset) ;
     % sigmoid activation function
    activation_teacher  = erf(teacher_scalar/sqrt(2));
    activation_student  = erf(student_scalar(j)/sqrt(2));
    % summation of all the dot products of student & teacher weight vectors
    % with input data.
    si = si + activation_student;
    tou = tou + activation_teacher;
    
    end
    % both layers getting updated simataneously for mu times.
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
       % updating the velocity using the gradient
        velocities = beta * velocities + (1-beta) * gradient_epsilon_studentweights';
         % updating the studen weight vector using updated velocity
        student_weights(j,:) =  student_weights(j,:) - learning_rate/samples*(velocities);
        
    end
    
    si = 0;
    tou = 0;
    
end


x = 0:mu/samples-1;
figure(1)

plot(x,R_in(:,1));
hold on;
plot(x,R_in(:,2));
plot(x,Q_ik(:,1));
plot(x,Q_ik(:,2));
xlabel('alpha')
ylabel('Order parameter momentum')
legend('R_{11}','R_{12}','Q_{11}','Q_{12}')
gen_error
figure(2)
plot(x, gen_error,'DisplayName', 'Momentum')
xlabel('alpha')
ylabel('Generalization error momentum')
legend('-DynamicLegend');
hold all