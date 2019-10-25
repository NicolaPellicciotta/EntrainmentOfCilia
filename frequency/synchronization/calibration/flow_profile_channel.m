%%%%% flow profile in channel with solution from White 1974
%%%%% channel squared 1 mm everything in mm
h=1.;
L=12.;
w=1.;
y=-w/2:0.01:w/2;
z= 0:0.01:h;
dp=1;
mu=1;
n_array=(1:2:20);
sumF=zeros([numel(z),numel(y)]);
for i=1:numel(n_array)
    n=n_array(i);
    temp= 1/(n^3)*(1- (cosh(n*pi*y./h)/cosh(n*pi*w/(2*h)))).*(sin(n*pi*z./h)');
    
    sumF=  sumF+ temp;
    
end

u= (4*(h^2)*dp)/((pi^3)*mu*L).* sumF;  

figure();
plot(y,u(2,:));hold on;
plot(y,u(floor(numel(z)/2),:));hold on;
xlabel('y');
ylabel('flow velocity')

 %% without dp and with Q flow rate
 
 %%% better
h=1;  %height of channel
w=1.; %width
y=-w/2:0.01:w/2;  % axis
z= 0:0.01:h;  %% height variable
Q= 50/60;   %% 10ul/min   mm3/s flow rate
n_even=(1:2:10);
n_disp=(2:2:10)-1;
 
sumF=zeros([numel(z),numel(y)]);
for i=1:numel(n_even)
    n=n_even(i);
    temp= 1/(n^3)*(1- (cosh(n*pi*y./h)/cosh(n*pi*w/(2*h)))).*(sin(n*pi*z./h)');
    
    sumF=  sumF+ temp;
    
end
 
  
 sumF2=0;
for i=1:numel(n_disp)
    n=n_disp(i);
    temp= (1- (1/(n^5))*(192/(pi^5))*(h/w)*tanh(n*pi/(2*h/w)));    

    sumF2=  sumF2+temp;
end


u= ((48*Q)/(pi*sumF2*w*h))*sumF; 

figure();%plot(y,u(floor(numel(z)/2),:));hold on;
plot(y,u(z==0.02,:));hold on;
xlabel('channel width');ylabel('flow velocity [mm/s]')

Q_exp=  sum(u(:))*0.01*0.01;
V_av= sum(u(:))/numel(u(:));