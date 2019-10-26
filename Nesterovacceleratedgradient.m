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
rngSeed = 999;
beta = 0.9 % maybe this scales 
%beta = 0.7666668/samples
epsilon = 2*pi
K=2
rng(5782)
teacher_weights = randn(K,samples)./sqrt(samples);
teacher_weights_ortho = GramSchmidt(teacher_weights');
teacher_weights = teacher_weights_ortho';
student_weights = randn(K,samples)./sqrt(samples);


velocities = 0;

for i = 1:mu
    % same array of input for the different weight vectors.
    % input to be generated for all the runs?
    
    
    if  learning_rate > 1.5
      %  learning_rate = learning_rate - 1/(samples*10)
       learning_rate = learning_rate - 0.1/samples
    end
    
   % if  beta > 0
   %   beta = cos(epsilon)*0.9;
   %   i;
   %   epsilon = epsilon - pi/10000;
   % else
   %     beta = 0;
   % end
    
    rng(rngSeed)
    dataset = randn(samples,1);
    rngSeed = rngSeed + 1;
    for j = 1:K
   %i
    
   % teacher_1 = dot(teacher_w1, dataset);
    %student_1 = dot(student_w1, dataset);
    teacher_scalar = dot(teacher_weights(j,:), dataset);
    %vW(t+1) = momentum.*Vw(t) - scaling .* gradient_F( W(t) + momentum.*vW(t) )
    student_scalar(j) = dot(student_weights(j,:)- beta* velocities, dataset) ;


    activation_teacher  = erf(teacher_scalar/sqrt(2));
    activation_student  = erf(student_scalar(j)/sqrt(2));
    
    si = si + activation_student;
    tou = tou + activation_teacher;
    
    end
    
    for j = 1:K
        if mod(i,samples) == 0
            R_in(i/samples,j) = dot(student_weights(1,:), teacher_weights(j,:));
            Q_ik(i/samples,j) = dot(student_weights(1,:), student_weights(j,:));
            gen_error(i/samples) = generalizationerror(student_weights,teacher_weights,K);
    end
        student_s = student_scalar(j);
        % gradient derivative done
        gradeint_epsilon_studentweights = (si - tou) * sqrt(2/pi)*exp(-(student_s*student_s)/2)* dataset;
       % velocities = beta * velocities + (1-beta)* ((learning_rate/samples).* gradeint_epsilon_studentweights');
        velocities = beta * velocities + (1-beta)* gradeint_epsilon_studentweights';
        student_weights(j,:) =  student_weights(j,:) - learning_rate/samples* (velocities);
        
    end
    
    si = 0;
    tou = 0;
    
end

%Q_ik = mat2gray(Q_ik);
%R_in = mat2gray(R_in);

x = 0:mu/samples-1;
figure(1)

plot(x,R_in(:,1));
hold on

plot(x,R_in(:,2));
%plot(x,R_in(:,3));
plot(x,Q_ik(:,1));
plot(x,Q_ik(:,2));
xlabel('alpha')
ylabel('Order parameter (Nesterov accelerated gradient)')
%plot(x,Q_ik(:,3));
%legend('R','R','R','Q','Q','Q')
legend('R','R','Q','Q')
gen_error
figure(2)
plot(x, gen_error,'m')
xlabel('alpha')
ylabel('Generalization error ')
legend('Nesterov accelerated gradient','Color',[0 0.9 0],'AutoUpdate','off')
hold on
