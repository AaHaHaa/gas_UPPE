j=6;

J = 0:20;
for i =1:length(J)
    cg(i) = ClebschGordan(j,J(i),2,0,0,0);
end
%cg(J==j) = 0;

figure;
plot(J,cg)