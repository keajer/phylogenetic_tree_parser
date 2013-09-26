%%% log likelihood
syms p1 p2 p3 n1 n2 n3 n4 n5

%likelihood function
l = (1/4*((p1*p2*p3)+(1-3*p1)*p2*p3 + p1*(1-3*p2)*p3 + p1*p2*(1-3*p3)))^n1 * ...
	(1/4*(1-3*p1)*(1-3*p2)*(1-3*p3) + 3/4*p1*p2*p3)^n2 *...
	(1/4*((1-3*p1)*p2*p3+p1*(1-3*p2)*(1-3*p3)) + 2/4*p1*p2*p3)^n3 *...
	(1/4*((1-3*p1)*(1-3*p2)*p3+p1*p2*(1-3*p3)) + 2/4*p1*p2*p3)^n4 *...
	(1/4*((1-3*p1)*p2*(1-3*p3)+p1*(1-3*p2)*p3) + 2/4*p1*p2*p3)^n5
	
%gradient

dp1 = diff(l, p1);

%log likelihood function
l1 = n1 * log(1/4*((p1*p2*p3)+(1-3*p1)*p2*p3 + p1*(1-3*p2)*p3 + p1*p2*(1-3*p3)));
l2 = n2 * log(1/4*(1-3*p1)*(1-3*p2)*(1-3*p3) + 3/4*p1*p2*p3);
l3 = n3 * log(1/4*((1-3*p1)*p2*p3+p1*(1-3*p2)*(1-3*p3)) + 2/4*p1*p2*p3);
l4 = n4 * log(1/4*((1-3*p1)*(1-3*p2)*p3+p1*p2*(1-3*p3)) + 2/4*p1*p2*p3);
l5 = n5 * log(1/4*((1-3*p1)*p2*(1-3*p3)+p1*(1-3*p2)*p3) + 2/4*p1*p2*p3);
	
dldp1 = diff(l1, p1) + diff(l2, p1) + diff(l3, p1) + diff(l4, p1) + diff(l5, p1);

%dldp1 = (n1*((11*p2*p3)/4 + p2*(3*p3 - 1) + p3*(3*p2 - 1))) / 
%(p1*p2*(3*p3 - 1) + p1*p3*(3*p2 - 1) + p2*p3*(3*p1 - 1) - (p1*p2*p3)/4)
% - (n3*((p2*p3)/4 - ((3*p2 - 1)*(3*p3 - 1))/4)) / 
%((p1*(3*p2 - 1)*(3*p3 - 1))/4 - (p2*p3*(3*p1 - 1))/4 + (p1*p2*p3)/2)
% + (n4*((p2*p3)/2 - p2*(3*p3 - 1) + (3*p3*(3*p2 - 1))/4))
%/((p3*(3*p1 - 1)*(3*p2 - 1))/4 - p1*p2*(3*p3 - 1) + (p1*p2*p3)/2)
% + (n5*((p2*p3)/2 + (3*p2*(3*p3 - 1))/4 - (p3*(3*p2 - 1))/4))
%/((p2*(3*p1 - 1)*(3*p3 - 1))/4 - (p1*p3*(3*p2 - 1))/4 + (p1*p2*p3)/2)
% + (n2*((3*p2*p3)/4 - (3*(3*p2 - 1)*(3*p3 - 1))/4))
% / ((3*p1*p2*p3)/4 - ((3*p1)/4 - 1/4)*(3*p2 - 1)*(3*p3 - 1))

dldp2 = diff(l1, p2) + diff(l2, p2) + diff(l3, p2) + diff(l4, p2) + diff(l5, p2);

%(n1*((11*p1*p3)/4 + p1*(3*p3 - 1) + p3*(3*p1 - 1)))/(p1*p2*(3*p3 - 1) + p1*p3*(3*p2 - 1) + p2*p3*(3*p1 - 1) - (p1*p2*p3)/4) - (n5*((p1*p3)/4 - ((3*p1 - 1)*(3*p3 - 1))/4))/((p2*(3*p1 - 1)*(3*p3 - 1))/4 - (p1*p3*(3*p2 - 1))/4 + (p1*p2*p3)/2) + (n3*((p1*p3)/2 + (3*p1*(3*p3 - 1))/4 - (p3*(3*p1 - 1))/4))/((p1*(3*p2 - 1)*(3*p3 - 1))/4 - (p2*p3*(3*p1 - 1))/4 + (p1*p2*p3)/2) + (n4*((p1*p3)/2 - p1*(3*p3 - 1) + (3*p3*(3*p1 - 1))/4))/((p3*(3*p1 - 1)*(3*p2 - 1))/4 - p1*p2*(3*p3 - 1) + (p1*p2*p3)/2) + (n2*((3*p1*p3)/4 - 3*((3*p1)/4 - 1/4)*(3*p3 - 1)))/((3*p1*p2*p3)/4 - ((3*p1)/4 - 1/4)*(3*p2 - 1)*(3*p3 - 1))

dldp3 = diff(l1, p3) + diff(l2, p3) + diff(l3, p3) + diff(l4, p3) + diff(l5, p3);

%(n1*((11*p1*p2)/4 + p1*(3*p2 - 1) + p2*(3*p1 - 1)))/(p1*p2*(3*p3 - 1) + p1*p3*(3*p2 - 1) + p2*p3*(3*p1 - 1) - (p1*p2*p3)/4) - (n4*((5*p1*p2)/2 - ((3*p1 - 1)*(3*p2 - 1))/4))/((p3*(3*p1 - 1)*(3*p2 - 1))/4 - p1*p2*(3*p3 - 1) + (p1*p2*p3)/2) + (n3*((p1*p2)/2 + (3*p1*(3*p2 - 1))/4 - (p2*(3*p1 - 1))/4))/((p1*(3*p2 - 1)*(3*p3 - 1))/4 - (p2*p3*(3*p1 - 1))/4 + (p1*p2*p3)/2) + (n5*((p1*p2)/2 - (p1*(3*p2 - 1))/4 + (3*p2*(3*p1 - 1))/4))/((p2*(3*p1 - 1)*(3*p3 - 1))/4 - (p1*p3*(3*p2 - 1))/4 + (p1*p2*p3)/2) + (n2*((3*p1*p2)/4 - 3*((3*p1)/4 - 1/4)*(3*p2 - 1)))/((3*p1*p2*p3)/4 - ((3*p1)/4 - 1/4)*(3*p2 - 1)*(3*p3 - 1))

p1hat = solve('dldp1', 'p1')

[p1hat, p2hat, p3hat] = solve('dldp1', 'dldp2', 'dldp3', 'p1', 'p2', 'p3')


