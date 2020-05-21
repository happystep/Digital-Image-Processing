% We will use this file to get used to the mechanics involved in matlab.
% Matlab will read the .m file as a script from the top.

% Variables are declared as follows don't forget statements end with ;

x = 3; % here we are assign the variable x with the value 
disp(x)

% Variables can be overwritten

x = 6; 
disp(x)


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

v = [1 4 7 10 13];

disp(v)

















 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 