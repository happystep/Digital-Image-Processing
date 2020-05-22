% We will use this file to get used to the mechanics involved in matlab.
% Matlab will read the .m file as a script from the top.

% Variables are declared as follows don't forget statements end with ;
% 1 - Tutorial Lessons 1

x = 3; % here we are assign the variable x with the value 
disp(x)

% Variables can be overwritten

x = 6; 
disp(x)

% 2 - Tutorial Lessons 2
% Matlab includes many bult in mathematical functions 
% elementary functions include, but are not limited to

disp(cos(x))
disp(sin(x))
disp(tan(x))
disp(abs(x))
disp(ceil(x))
disp(floor(x))

% Matlab also includes predefined constant values

disp(pi)
disp(1i)
disp(Inf)
disp(NaN)

% Example using defined elementary functions  previously defined.

 a = 5; x = 2; y = 8;
 y = exp(-a)*sin(x)+10*sqrt(y);
 disp(y)
 
 disp(log(142))
 disp(log10(142))
 
 disp(sin(pi/4))
 disp(exp(10))
 
 % next we will dive into plotting using the command plot(x,y)
 
 x = [1 2 3 4 5 6];
 y = [3 -1 2 4 5 1];
 plot(x,y)
 
% For example, to plot the function sin (x) on the interval [0, 2π], we first create a vector of
% x values ranging from 0 to 2π, then compute the sine of these values, and finally plot the result:

x = 0:pi/100:2*pi;
y = sin(x);
plot(x,y)

% we can add titles, labels and annotations to the figure

xlabel('x = 0:2\pi')
ylabel('Sine of x')
title('Plot of the Sine function')
 

% multiple data sets in one plot

x = 0:pi/100:2*pi;
y1 = 2*cos(x);
y2 = cos(x);
y3 = 0.5*cos(x);
plot(x,y1,'--',x,y2,'-',x,y3,':')
xlabel('0 \leq x \leq 2\pi')
ylabel('Cosine functions')
legend('2*cos(x)','cos(x)','0.5*cos(x)')
title('Typical example of multiple plots')
axis([0 2*pi -3 3])

% specifying line styles and colors

% plot(x,y,'style_color_marker') % where style color marker is a triplet of
% color, symbol line style and symbol marker eg: k-+ or m:s

% Matrix generalization, you have row vectors or column vectors. Matrices
% are a two-dimensional array consisting of m rows and n columns. 

% row vector
v = [1 4 7 10 13];
disp(v)

% Column vectors are created in a similar way, however, semicolon (;) 
% must separate the components of a column vector

w = [1;4;7;10;13];
disp(w)

% one can also transpose a row vector to convert it into a column vector

w = v';
disp(w)

% v(1) is the first element of vector v, v(2) its second element, and so forth.
%Furthermore, to access blocks of elements, we use MATLAB’s colon notation (:). 
% For example, to access the first three elements of v, we write

disp(v(1:3))

disp(v(:)) % produces a column vector

disp(v(1:end)) % produces a row vector

% entering and displaying a matrix

A = [1 2 3; 4 5 6; 7 8 9];

disp(A)
disp (A(2,1)) % (row,col)

% matrix indexing can be used to replace a value of an existing matrix

A(3,3) = 0;
disp(A)

% The colon operator will prove very useful and understanding how it works is the key to
% efficient and convenient usage of MATLAB. It occurs in several different forms.
% Often we must deal with matrices or vectors that are too large to enter one element at a time.
% For example, suppose we want to enter a vector x consisting of points
% (0, 0.1, 0.2, 0.3, · · · , 5). We can use the command

x = 0:0.1:5;

disp(x)

% linear spacing

theta = linspace(0, 2*pi, 101);
disp(theta)

% colon operator in a matrix

disp(A(2,:)) % second row of elements in A

disp(A(:,2:3)) % sub matrix with the last two columns of A

% A row or a column of a matrix can be deleted by setting it to a null vector, [ ].

%A(:,2)=[];
%disp(A)

% creating a sub-matrix

% To extract a submatrix B consisting of rows 2 and 3 and columns 1 and 2 of the matrix A, do the following
B = A([2 3],[1 2]);
disp(B)


% To interchange rows 1 and 2 of A, use the vector of row indices together with the colon operator
C = A([2 1 3],:);
disp(C)

% To create a vector version of matrix A, do the following

disp(A(:))

% deleting and restoring

 A(3,:) = [];
 A = [A(1,:);A(2,:);[7 8 0]];
 disp(A)
 
 % dimension
 
 disp(size(A))
 
 % continuation
 
 % If it is not possible to type the entire input on the same line, use consecutive periods, called
 % an ellipsis . . ., to signal continuation, then continue the input on the next line.
 
 % B = [4/5 7.23*tan(x) sqrt(6); ...
 %       1/x.^2 0 3/(x*log(x)); ...
 %      x-7 sqrt(3) x*sin(x)];
 
% disp(B)

% Transposing a matrix is the same as transposing a vector

disp(A')

% concatenating matrices, matrices made up of sub-matrices

B = [A 10*A; -A [1 0 0; 0 1 0; 0 0 1]];

disp(B)


% elementary matrices built in

n = 3;
m = 3;

disp(eye(m,n))
disp(eye(n))
disp(zeros(m,n)) 
disp(ones(m,n))
disp(diag(A))

% Two other important matrix generation functions are rand and randn, which generate
% matrices of (pseudo-)random numbers using the same syntax as eye.
disp(rand(m,n))

% In addition, matrices can be constructed in a block form. With C defined by C = [1 2; 3 4], we may create a matrix D as follows
C = [1 2; 3 4];
D = [C zeros(2); ones(2) eye(2)];
disp(D)


% Special Matrices 

% hilb Hilbert matrix
% invhilb Inverse Hilbert matrix
% magic Magic square
% pascal Pascal matrix
% toeplitz Toeplitz matrix
% vander Vandermonde matrix
% wilkinson Wilkinson’s eigenvalue test matrix


% 3 - Array operations and linear equations

% matrix arithmetic operations, +, -, *, ^ all work

% array arithmetic operations, + and - work the same as matrix


A = [1 2 3; 4 5 6; 7 8 9];
B = [10 20 30; 40 50 60; 70 80 90];

C = A.*B;
disp(C)

D = A.^2;
disp(D)
 
 
% solving linear equations

A = [1 2 3; 4 5 6; 7 8 0];
b = [1; 1; 1];
x = inv(A)*b;
disp(x)

% or this way, is faster computationally
A = [1 2 3; 4 5 6; 7 8 0];
b = [1; 1; 1];
x = A\b;
disp(x)

% we also have inverse and determinants
disp(inv(A))
disp(det(A))

% matrix functions include

disp(eig(A))
disp(norm(A))
disp(rank(A))

% 4 - Introduction to Programming with MATLAB

% f = factorial(5);
% disp(f)

% use of input
num = input('Enter a number between 0-10: ');
disp(num^2)

% function f = factorial(n) 
% % FACTORIAL(N) returns the factorial of N. 
% % Compute a factorial value. 
% f = prod(1:n);
% end

% function is commented out in order to continue with book

% 5 - Control flow and operators

% control flow - if, elseif, else, end
a = 1; 
b = 2;
c = 3;
discr = b*b - 4*a*c;
if discr < 0
disp('Warning: discriminant is negative, roots are imaginary')
end

discr = b*b - 4*a*c;
if discr < 0
disp('Warning: discriminant is negative, roots are imaginary');
else
disp('Roots are real, but may be repeated')
end


discr = b*b - 4*a*c;
if discr < 0
disp('Warning: discriminant is negative, roots are imaginary');
elseif discr == 0
disp('Discriminant is zero, roots are repeated')
else
disp('Roots are real')
end

% control flow - rational and logical operators

% > Greater than
% < Less than
% >= Greater than or equal to
% <= Less than or equal to
% == Equal to
% ∼= Not equal to
% & AND operator
% | OR operator
% ∼ NOT operator


% control flow - for end loop

for ii=1:5
x=ii*ii;
disp(x)
end

n = 5; A = eye(n);
for j=2:n
for i=1:j-1
A(i,j)=i/j;
A(j,i)=i/j;
end
end

disp(A)
 
% control flow - while end loop

x = 1;
while x <= 10
x = 3*x;
disp(x)
end


% other flow strucutre include, break, continue and return. 

% MATLAB Operator Precedence

% 1 Parentheses
% 2 Transpose 
% 3 Unary plus (+), unary minus (−), logical negation (∼)
% 4 Multiplication (. ∗), right division (. /), left division (.\),
%     matrix multiplication (∗), matrix right division (/),
%     matrix left division (\)
% 5 Addition (+), subtraction (−)
% 6 Colon operator (:)
% 7 Less than (<), less than or equal to (≤), greater (>),
% greater than or equal to (≥), equal to (==), not equal to (∼=)
% 8 Element-wise AND, (&)
% 9 Element-wise OR, (|)



% writing to a file

% write some variable length strings to a file
op = fopen('weekdays.txt','wt');
fprintf(op,'Sunday\nMonday\nTuesday\nWednesday\n');
fprintf(op,'Thursday\nFriday\nSaturday\n');
fclose(op);
 
% 6 - Debugging M-files

% this is similar to other IDE. The user can set different breakpoints. 
 
 
 