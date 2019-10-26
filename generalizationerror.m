function [generalization_error] = generalizationerror(student_weights,teacher_weights,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
student_weights = student_weights;
teacher_weights = teacher_weights;
k=k;
sum1 = 0;
sum2 = 0;
sum3 = 0;

for i = 1:k
    for j = 1:k
        sum1 = sum1 + asin((dot(student_weights(i,:),student_weights(j,:)))/(sqrt(1+dot(student_weights(i,:),student_weights(i,:)))*sqrt(1+dot(student_weights(j,:),student_weights(j,:)))));
        sum2 = sum2+  asin((dot(teacher_weights(i,:),teacher_weights(j,:)))/(sqrt(1+dot(teacher_weights(i,:),teacher_weights(i,:)))*sqrt(1+dot(teacher_weights(j,:),teacher_weights(j,:)))));
        sum3 = sum3 + asin((dot(student_weights(i,:),teacher_weights(j,:)))/(sqrt(1+dot(student_weights(i,:),student_weights(i,:)))*sqrt(1+dot(teacher_weights(j,:),teacher_weights(j,:)))));
    end
end

generalization_error = (1/pi) * (sum1+sum2-2*sum3);

end

