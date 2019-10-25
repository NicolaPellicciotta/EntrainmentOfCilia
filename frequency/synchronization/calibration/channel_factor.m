function F = channel_factor(z, y, w, h)
% INPUT: a z value for which to calculate the factor, y position, width and height of
% the channel
%OUTPUT: the numerical factor that depends on channel geometry, from
%White's solution

term1 = sin(pi*z/h) * (1 - cosh(pi*y/h)/cosh(pi*w/(2*h)));
term2 = (1/27) * sin(3*pi*z/h) * (1 - cosh(3*pi*y/h)/cosh(3*pi*w/(2*h)));
term3 = (192*h/(pi^5*w))*tanh(pi*w/2*h);
term4 = (192*h/(243*pi^5*w))*tanh(3*pi*w/(2*h));
factor1 = 48/(pi^3*h*w);

F = factor1 * (term1 + term2)/(1 - (term3 + term4));

end