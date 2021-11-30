
function C=C_analyse(M,q1,q2,dq2,dq1)
%Christoffel

c111=simplify(1/2*(diff(M(1,1),q1)+diff(M(1,1),q1)-diff(M(1,1),q1)));
c112=simplify(1/2*(diff(M(2,1),q1)+diff(M(2,1),q1)-diff(M(1,1),q2)));
c121=simplify(1/2*(diff(M(1,2),q1)+diff(M(1,1),q2)-diff(M(1,2),q1)));
c122=simplify(1/2*(diff(M(2,2),q1)+diff(M(2,1),q2)-diff(M(1,2),q2)));
c211=c121;                                                           
c212=c122;                                                           
c221=simplify(1/2*(diff(M(1,2),q2)+diff(M(1,2),q2)-diff(M(2,2),q1)));
c222=simplify(1/2*(diff(M(2,2),q2)+diff(M(2,2),q2)-diff(M(2,2),q2)));

C=[c121*dq2  c211*dq2+c221*dq1
   -c112*dq1            0];

end