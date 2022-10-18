sig = @(v)(4*(exp(2.5*(v-1))-1));
v = 0:0.001:1.2;
plot(v,sig(v))

dsig = @(E)(10+2.5*E);
E = 0:50;
plot(E,dsig(E))

lnsig = @(v)(log(4*(exp(2.5*(v-1))-1)));
plot(v,lnsig(v))
plot(log(v), sig(v))