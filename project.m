close all
clc;clear;
rng default
 
%create w1 samples
point_a = 2; point_b = 8;
point_c = 1; point_d= 2;
x1 = point_a + (point_b-point_a).*rand(400,1);
x2 = point_c + (point_d-point_c).*rand(400,1);
 
w1 = [x1 ,x2];
 
%create w2 samples
point_a = 6; point_b = 8;
point_c = 2.5; point_d= 5.5;
x1 = point_a + (point_b-point_a).*rand(100,1);
x2 = point_c + (point_d-point_c).*rand(100,1);
 
w2 = [x1 ,x2];
 
%plot samples
figure(1)
plot(w1(:,1),w1(:,2),'b.', w2(:,1),w2(:,2),'gx', 0,0);
title('2-D Sample Plot');
legend('w1 Sample','w2 Sample','Location','northwest');
 


%Meros B
%B1
%find mean values of samples
m1 = mean(w1);
m2 = mean(w2);

%find covariances
s1 = cov(w1);
s2 = cov(w2);
 
p11 = [];
p12 = [];
for i = 1:400
    p11 = [p11,(1/(2*pi*sqrt(abs(det(s1))))) * exp(-(1/2)*(w1(i,:) - m1)*inv(s1)*(w1(i,:) - m1)')];
    p12 = [p12,(1/(2*pi*sqrt(abs(det(s2))))) * exp(-(1/2)*(w1(i,:) - m2)*inv(s2)*(w1(i,:) - m2)')];
    
end
 
p21 = [];
p22 = [];
for i = 1:100
    p22 = [p22,(1/(2*pi*sqrt(abs(det(s2))))) * exp(-(1/2)*(w2(i,:) - m2)*inv(s2)*(w2(i,:) - m2)')];
    p21 = [p21,(1/(2*pi*sqrt(abs(det(s1))))) * exp(-(1/2)*(w2(i,:) - m1)*inv(s1)*(w2(i,:) - m1)')];
end


%B2
w1_right = [];
w2_right = [];
w1_wrong = [];
w2_wrong = [];

%Taksinomisi stoixeiwn W1 me euclidean distance
for i = 1:400
    d1 = norm(w1(i,:) - m1);
    d2 = norm(w1(i,:) - m2);
    if d1<=d2
        w1_right = [w1_right;w1(i,:)];
    else
        w1_wrong = [w1_wrong;w1(i,:)];
    end
end

%Taksinomisi stoixeiwn W2 me euclidean distance
for i = 1:100
    d1 = norm(w2(i,:) - m1);
    d2 = norm(w2(i,:) - m2);
    if d1<d2
        w2_wrong = [w2_wrong;w2(i,:)];
    else
        w2_right = [w2_right;w2(i,:)];
    end
end
 
%Ypologismos sfalmatos taksinomisis
error = (size(w1_wrong,1) + size(w2_wrong,1)) / 500;
disp(['Euclidean distance error : ',num2str(error*100),'%']);
 
%Plot dedomena klasewn & lathi euclidean taksinomisis
figure(2)
plot(w1(:,1),w1(:,2),'b.',w2(:,1),w2(:,2),'gx', 0,0);
 
if size(w1_wrong,1) ~= 0 
    hold on;
    plot(w1_wrong(:,1),w1_wrong(:,2),'r.', 0,0);
end

if size(w2_wrong,1) ~= 0
    hold on;
    plot(w2_wrong(:,1),w2_wrong(:,2),'rx', 0,0);
end

% Plot mesa dianysmata klasewn
plot(m1(:,1), m1(:,2), 'k^', 0,0);
hold on;
plot(m2(:,1), m2(:,2), 'k^', 0,0);

 
%B3
s = (0.8*s1)+(0.2*s2);
w1_right = [];
w2_right = [];
w1_wrong = [];
w2_wrong = [];

%Taksinomisi stoixeiwn W1 me mahalanobis distance 
for i = 1:400
    d1 = sqrt((w1(i,:) - m1)*inv(s)*(w1(i,:) - m1)');
    d2 = sqrt((w1(i,:) - m2)*inv(s)*(w1(i,:) - m2)');
    if d1<=d2
        w1_right = [w1_right;w1(i,:)];
    else
        w1_wrong = [w1_wrong;w1(i,:)];
    end
end

%Taksinomisi stoixeiwn W2 me mahalanobis distance
for i = 1:100
    d1 = sqrt((w2(i,:) - m1)*inv(s)*(w2(i,:) - m1)');
    d2 = sqrt((w2(i,:) - m2)*inv(s)*(w2(i,:) - m2)');
    if d1<d2
        w2_wrong = [w2_wrong;w2(i,:)];
    else
        w2_right = [w2_right;w2(i,:)];
    end
end
 
%Ypologismos sfalmatos mahalanobis distance
error = (size(w1_wrong,1) + size(w2_wrong,1)) / 500;
disp(['Mahalanobis distance error : ',num2str(error*100),'%']);

%Plot lathi mahalanobis distance
if size(w1_wrong,1) ~= 0 
    %hold on;
    plot(w1_wrong(:,1),w1_wrong(:,2),'k.', 0,0);
end

if size(w2_wrong,1) ~= 0
    %hold on;
    plot(w2_wrong(:,1),w2_wrong(:,2),'kx', 0,0);
end

title('2-D Sample Plot with errors');
legend('w1 Sample','w2 Sample','Euc. errors w1','Euc. errors w2', 'Mah. erros w1','Mah. erros w2','Location','northwest');
 

%B4
w1_right = [];
w2_right = [];
w1_wrong = [];
w2_wrong = [];

for i = 1:400
    if p11>=p12
        w1_right = [w1_right;w1(i,:)];
    else
        w1_wrong = [w1_wrong;w1(i,:)];
    end
end

for i = 1:100
    if p21>p22
        w2_wrong = [w2_wrong;w2(i,:)];
    else
        w2_right = [w2_right;w2(i,:)];
    end
end
 
error = (size(w1_wrong,1) + size(w2_wrong,1)) / 500;
disp(['Bayesian classification error : ',num2str(error*100),'%']);
 


% Meros C
%C1
figure(3);

projectedPCAw1 = w1;
projectedPCAw2 = w2;

%Pinakas Dedomenwn
M = [w1;w2];

%Covariance Pinaka
Matrix = cov(M);

%Idiotimes
eigen = eig(Matrix);

%Max idiotimi
eigen_max = max(eigen);

colV = [1;1];
colV = eigen_max*colV;

%Evresi eytheias diaxwrismou
z = Matrix\colV;

syms x y;

% Plot eytheias diaxwrismou & deigmatwn
f1(x,y) = z(1,:)*x + z(2,:)*y;
h = ezplot(f1,[-4,10]);
set(h, 'Color', 'k');
hold on
plot(w1(:,1),w1(:,2),'b.');
hold on;
plot(w2(:,1),w2(:,2),'gx');

title('PCA');

% Plot eytheias diaxwrismou
figure
f1(x,y) = z(1,:)*x + z(2,:)*y;
h = ezplot(f1,[-4,10]);
set(h, 'Color', 'k');
hold on

slope = -z(1,:)/z(2,:);

%Evresi projected points gia W1
for i = 1:400
   b0 = -(-1/slope)*w1(i,1) + w1(i,2);
   xIntersection = (b0)/ (slope - (-1/slope));
   yIntersection = (-1/slope) * xIntersection + b0;

   projectedPCAw1(i,1)=xIntersection;
   projectedPCAw1(i,2)=yIntersection;
    
   plot(xIntersection, yIntersection, 'b.')

end

%Evresi projected points gia W2
for i = 1:100
   b0 = -(-1/slope)*w2(i,1) + w2(i,2);
   xIntersection = (b0)/ (slope - (-1/slope));
   yIntersection = (-1/slope) * xIntersection + b0;
    
   projectedPCAw2(i,1)=xIntersection;
   projectedPCAw2(i,2)=yIntersection;
    
   plot(xIntersection, yIntersection, 'gx')

end
title('PCA - PROJECTIONS')



%C2
w1_wrong = [];
w1_right = [];
w2_wrong = [];
w2_right = [];

mesos_projectedPCA_w1 = mean(projectedPCAw1);
mesos_projectedPCA_w2 = mean(projectedPCAw2);

figure

%Plot eytheias diaxwrismou
f1(x,y) = z(1,:)*x + z(2,:)*y;
h = ezplot(f1,[-4,10]);
set(h, 'Color', 'k');
hold on

%Evresi projected points gia W1
for i = 1:400
    
   b0 = -(-1/slope)*w1(i,1) + w1(i,2);
   xIntersection = (b0)/ (slope - (-1/slope));
   yIntersection = (-1/slope) * xIntersection + b0;

   projectedPCAw1(i,1)=xIntersection;
   projectedPCAw1(i,2)=yIntersection;
    
   plot(xIntersection, yIntersection, 'b.') 
    
   %Taksinomisi projected points
   d1 = norm(projectedPCAw1(i,:) - mesos_projectedPCA_w1);
   d2 = norm(projectedPCAw1(i,:) - mesos_projectedPCA_w2);
   if d1<=d2
        w1_right = [w1_right;projectedPCAw1(i,:)];
   else
        w1_wrong = [w1_wrong;projectedPCAw1(i,:)];
   end
end

%Evresi projected points gia W2
for i = 1:100

   b0 = -(-1/slope)*w2(i,1) + w2(i,2);
   xIntersection = (b0)/ (slope - (-1/slope));
   yIntersection = (-1/slope) * xIntersection + b0;
    
   projectedPCAw2(i,1)=xIntersection;
   projectedPCAw2(i,2)=yIntersection;
    
   plot(xIntersection, yIntersection, 'gx')
    
   %Taksinomisi projected points
   d1 = norm(projectedPCAw2(i,:) - mesos_projectedPCA_w1);
   d2 = norm(projectedPCAw2(i,:) - mesos_projectedPCA_w2);
   if d1<d2
       w2_wrong = [w2_wrong;projectedPCAw2(i,:)];
   else
       w2_right = [w2_right;projectedPCAw2(i,:)];
   end
end

%Plot PCA me lathi taksinomisis
title('PCA errors');
if size(w1_wrong,1) ~= 0 
    hold on;
    plot(w1_wrong(:,1),w1_wrong(:,2),'r.', 0,0);
end
if size(w2_wrong,1) ~= 0
    hold on;
    plot(w2_wrong(:,1),w2_wrong(:,2),'rx', 0,0);
end

disp(['Pososto lathous PCA : ',num2str( ((size(w1_wrong,1)+size(w2_wrong,1))/500)*100 ),'%']);


%C3
figure

% P(W1) & P(W2)
Pithanotita_W1 = 400/500;
Pithanotita_W2 = 100/500;

% Pinakes dedomenwn
projectedLDAw1 = w1;
projectedLDAw2 = w2;

Pinakas_1_lda = w1-m1;
Pinakas_2_lda = w2-m2;

Pinakas_1_lda = Pithanotita_W1*(transpose(Pinakas_1_lda)*Pinakas_1_lda);
Pinakas_2_lda = Pithanotita_W2*(transpose(Pinakas_2_lda)*Pinakas_2_lda);

Sw = Pithanotita_W1*Pinakas_1_lda+Pithanotita_W2*Pinakas_2_lda;
wproj = inv(Sw)*(transpose(m1)-transpose(m2));

slope_lda = -wproj(1,:)/wproj(2,:);
slope2 = -wproj(1,:)/wproj(2,:);

syms x y;

%Plot LDA eytheias diaxwrismou & klasewn
f3(x,y) = (-1/slope2)*x + y;
h = ezplot(f3,[-4,10]);
set(h, 'Color', 'k');
hold on
plot(w1(:,1),w1(:,2),'b.');
hold on;
plot(w2(:,1),w2(:,2),'gx');

title('LDA');

%Plot eytheias diaxwrismou
figure
f3(x,y) = (-1/slope2)*x + y;
h = ezplot(f3,[-4,10]);
set(h, 'Color', 'k');
hold on

%Evresi projected points gia W1
for i = 1:400
   
   b0 = -(slope2)*w1(i,1) + w1(i,2);
   xIntersection = -(b0)/ ((slope2)-(-1/slope2));
   yIntersection = (slope2) * xIntersection + b0;

   projectedLDAw1(i,1)=xIntersection;
   projectedLDAw1(i,2)=yIntersection;
  
   plot(xIntersection, yIntersection, 'b.')

end

%Evresi projected points gia W2
for i = 1:100
   
   b0 = -(slope2)*w2(i,1) + w2(i,2);
   xIntersection = -(b0)/ ((slope2)-(-1/slope2));
   yIntersection = (slope2) * xIntersection + b0;

   projectedLDAw2(i,1)=xIntersection;
   projectedLDAw2(i,2)=yIntersection;
   
   plot(xIntersection, yIntersection, 'gx')
    
end
title('LDA - PROJECTIONS')

 
%C4
w1_wrong = [];
w1_right = [];
w2_wrong = [];
w2_right = [];

mesos_projectedLDA_w1 = mean(projectedLDAw1);
mesos_projectedLDA_w2 = mean(projectedLDAw2);

figure

% Plot eytheia diaxwrismou
f3(x,y) = -(1/slope2)*x + y;
h = ezplot(f3,[-4,10]);
set(h, 'Color', 'k');
hold on

% Evresi projected points W1
for i = 1:400
    
   b0 = -(slope2)*w1(i,1) + w1(i,2);
   xIntersection = -(b0)/ ((slope2)-(-1/slope2));
   yIntersection = (slope2) * xIntersection + b0;

   projectedLDAw1(i,1)=xIntersection;
   projectedLDAw1(i,2)=yIntersection;
    
   plot(xIntersection, yIntersection, 'b.') 
    
   % Taksinomisi projected points
   d1 = norm(projectedLDAw1(i,:) - mesos_projectedLDA_w1);
   d2 = norm(projectedLDAw1(i,:) - mesos_projectedLDA_w2);
   if d1<=d2
        w1_right = [w1_right;projectedLDAw1(i,:)];
   else
        w1_wrong = [w1_wrong;projectedLDAw1(i,:)];
   end
end

%Evresi projected points W2
for i = 1:100

   b0 = -(slope2)*w2(i,1) + w2(i,2);
   xIntersection = -(b0)/ ((slope2)- (-1/slope2));
   yIntersection = (slope2) * xIntersection + b0;
    
   projectedLDAw2(i,1)=xIntersection;
   projectedLDAw2(i,2)=yIntersection;
    
   plot(xIntersection, yIntersection, 'gx')
    
   %Taksinomisi projected points
   d1 = norm(projectedLDAw2(i,:) - mesos_projectedLDA_w1);
   d2 = norm(projectedLDAw2(i,:) - mesos_projectedLDA_w2);
   if d1<d2
       w2_wrong = [w2_wrong;projectedLDAw2(i,:)];
   else
       w2_right = [w2_right;projectedLDAw2(i,:)];
   end
end

%Plot me lathi taksinomisis
title('LDA errors');
if size(w1_wrong,1) ~= 0 
    hold on;
    plot(w1_wrong(:,1),w1_wrong(:,2),'r.', 0,0);
end
if size(w2_wrong,1) ~= 0
    hold on;
    plot(w2_wrong(:,1),w2_wrong(:,2),'rx', 0,0);
end

disp(['Pososto lathous LDA : ',num2str( ((size(w1_wrong,1)+size(w2_wrong,1))/500)*100 ),'%']);

% Meros D
%D1 Minimum Least Squares
y = ones(500,1);
y(401:500) = -y(401:500);
X = ones(500,3);
X(:,1:2) = [w1;w2];
w = inv(X'*X)*X'*y;
 
%Plot eytheias diaxwrismou & klasewn
figure;
plot(w1(:,1),w1(:,2),'b.',w2(:,1),w2(:,2),'gx', 0,0);
hold on;
syms x y;
f(x,y) = w(1)*x + w(2)*y + w(3);
ezplot(f,[0,10]);
title('Least Squares');
 
%Evresi sfalmatos tetragwnwn
J = [];
y = ones(500,1);
for i = 1:500
    J = [J, (y(i) - X(i,:)*w)^2];
end
lse = sum(J);
disp(['Minimum least Square, Sfalma Tetragwnwn: ',num2str(lse)]);

%Evresi sfalmatos taksinomisis
errors = 0;
for i = 1:400
    if w'*[w1(i,:),1]'<=0
        errors = errors +1;
    end
end
for  i = 1:100
    if w'*[w2(i,:),1]' >= 0 
        errors = errors +1;
    end
end

disp(['Minimum least Square, Error: ',num2str( (errors/500)*100 ),'%']);
 

%D2 Perceptron 
w = [1; 0; 7];
r = 1/50;
 
y = ones(500,1);
y(401:500) = -y(401:500);
X = [w1;w2];
X = [X,ones(500,1)];
 
Y = []; 
e = 0;
a = 1;

syms x y;

%Epanaliptiko loop taksinomisi twn dedomenwn
%kai epanaprosdiorismou toy w
while length(Y)|| a ~= 0
    a = 0;
    Y = [];
    for i = 1:500
        if i<=400
            if X(i,:)*w<=0
                Y = [Y;1*X(i,:)];
            end
        else
            if X(i,:)*w>=0
                Y = [Y;(-1)*X(i,:)];
            end
        end
    end
    e = e +1;
    sm = [];
    for i = 1:size(Y,1)
      sm = [sm ;Y(i,:)];
    end
    a = sum(sm);
    w = w +r*a';
 
end

%Plot eytheias diaxwrismou & klasewn
f(x,y) = w(1)*x + w(2)*y+w(3);
figure
plot(w1(:,1),w1(:,2),'b.',w2(:,1),w2(:,2),'gx', 0,0);
hold on;
xx = linspace(0,10,1000);
plot(xx,(xx*w(1)+w(3))/(-w(2)), 0,0)
legend('w1','w2','g(x)','Location', 'northwest');
title('Perceptron');
disp(['Number of Iterations : ',num2str(e)]);
