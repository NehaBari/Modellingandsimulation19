%% The script plots the orderparameter and generalization error for various learning rates.
% The orderparameter plot can get clumsy but the aim of this script is to
% display different GE for various learning rates
student_scalar = [];
teacher_scalar = [];
gen_error = [];
R_in = [];
Q_ik = [];
mu = 140000;
si = 0;
tou = 0;
samples= 1000;
K=2
rngSeed = 999;
Learningrates_array = [0.5,1,1.5,10];
sizeofleare=numel(Learningrates_array);
teacher_weights = randn(K,samples)./sqrt(samples); % teacher weight vector
teacher_weights_ortho = GramSchmidt(teacher_weights'); %conversion to Isotropic teachers
teacher_weights = teacher_weights_ortho';
student_weights = randn(K,samples)./sqrt(samples);
%velocities = randn(K,samples)./sqrt(samples);
velocities = 0;
for learning = 1:sizeofleare
    Learningrates_array(learning)
for i = 1:mu
    % same array of input for the different weight vectors.
    % input to be generated for all the runs?
    
    
    rng(rngSeed);
    dataset = randn(samples,1);
    rngSeed = rngSeed + 1;
    for j = 1:K

    teacher_scalar = dot(teacher_weights(j,:), dataset);
    student_scalar(j) = dot(student_weights(j,:), dataset) ;

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
        % updating the studen weight vector using gradient
        student_weights(j,:) =  student_weights(j,:) - ((Learningrates_array(learning)/samples)* gradient_epsilon_studentweights');

        
    end
    
    si = 0;
    tou = 0;
    
end
x = 0:mu/samples-1;
figure(1)
plot(x,R_in(:,1));
hold on
plot(x,R_in(:,2));
plot(x,Q_ik(:,1));
plot(x,Q_ik(:,2));
xlabel('alpha')
ylabel('Order-Parameters')
legend('R_{11}','R_{12}','Q_{11}','Q_{12}')
gen_error;
figure(2)
plot(x, gen_error)
xlabel('alpha')
ylabel('Generalization error')
legend('Learning rate optimization0.5','Learning rate optimization1','Learning rate optimization1.5','Learning rate optimization10')
hold on
end


