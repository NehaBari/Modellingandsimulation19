%%
student_scalar = [];
teacher_scalar = [];
gen_error = [];
R_in = [];
Q_ik = [];
mu = 140000;
si = 0;
tou = 0;
samples= 1000;
learning_rate = 1.5;
K=2
rngSeed = 999;
rng(5782)
teacher_weights = randn(K,samples)./sqrt(samples);
teacher_weights_ortho = GramSchmidt(teacher_weights');
teacher_weights = teacher_weights_ortho';
% x = ones(1,samples/2);
% y = zeros(1,samples/2);
% 
% x1 = [x,y]
% x2 = [y,x]
% 
% teacher_weights = [x1;x2]./sqrt(samples);
%dot(x1,x2)

student_weights = randn(K,samples)./sqrt(samples);
%velocities = randn(K,samples)./sqrt(samples);
velocities = 0;

for i = 1:mu
    % same array of input for the different weight vectors.
    % input to be generated for all the runs?
    
    
    if  learning_rate > 1.5
      %  learning_rate = learning_rate - 1/(samples*10)
       learning_rate = learning_rate - 0.1/(samples*0.5)
    end
    
    rng(rngSeed)
    dataset = randn(samples,1);
    rngSeed = rngSeed + 1;
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
    % both layers getting updated simataneously for mu times.
    for j = 1:K
        if mod(i,samples) == 0
            R_in(i/samples,j) = dot(student_weights(1,:), teacher_weights(j,:));
            Q_ik(i/samples,j) = dot(student_weights(1,:), student_weights(j,:));
            % gen_error = vector for number of steps. generalisation error
            % is for both the weight vectors together in one vector.
            gen_error(i/samples) = generalizationerror(student_weights,teacher_weights,K);
    end
        student_s = student_scalar(j);
        
        gradeint_epsilon_studentweights = (si - tou) * sqrt(2/pi)*exp(-(student_s*student_s)/2)* dataset;
       % student_weights(j,:) =  student_weights(j,:) - ((learning_rate/samples)* gradeint_epsilon_studentweights');
       student_weights(j,:) =  student_weights(j,:) - ((learning_rate/samples)* gradeint_epsilon_studentweights');
         %velocities = beta * velocities + (1-beta)* ((learning_rate/samples).* gradeint_epsilon_studentweights');
         %student_weights(j,:) =  student_weights(j,:) - (0.5*velocities);
        
    end
    
    si = 0;
    tou = 0;
    
end

%Q_ik = mat2gray(Q_ik);
%R_in = mat2gray(R_in);
% if wanna generate the order parameter graphs for one rate that is best
% then need to create a different script probably
x = 0:mu/samples-1;
figure(1)

plot(x,R_in(:,1));
hold on

plot(x,R_in(:,2));
%plot(x,R_in(:,3));
plot(x,Q_ik(:,1));
plot(x,Q_ik(:,2));
%plot(x,Q_ik(:,3));
%legend('R','R','R','Q','Q','Q')
xlabel('alpha')
ylabel('Order-Parameters-rate')
legend('R','R','Q','Q')
gen_error;
figure(2)
plot(x, gen_error)
xlabel('alpha')
ylabel('Generalization error rate')
%legend('LR 10')
hold on

