function fb = tmp(in1,in2,in3,in4,in5)
%TMP
%    FB = TMP(IN1,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    29-Mar-2024 21:54:43

n1 = in1(:,1);
n2 = in1(:,2);
param1 = in5(:,1);
u1 = in2(:,1);
u2 = in2(:,2);
u3 = in2(:,3);
u4 = in2(:,4);
uh1 = in3(:,1);
uh2 = in3(:,2);
uh3 = in3(:,3);
uh4 = in3(:,4);
uinf1 = in4(:,1);
uinf2 = in4(:,2);
uinf3 = in4(:,3);
uinf4 = in4(:,4);
t2 = uh2.^2;
t3 = uh3.^2;
t4 = 1.0./param1;
t5 = 1.0./uh1;
t7 = uh1.*uh4.*2.0;
t6 = t2+t3-t7;
t8 = 1.0./t6;
t9 = 1.0./uh1.^2;
t10 = param1-1.0;
t11 = sqrt(2.0);
t12 = n1.*uh2.*2.0;
t13 = n2.*uh3.*2.0;
t17 = param1.*t6.*t9.*t10;
t14 = sqrt(-t17);
t15 = t11.*t14.*uh1;
t16 = t12+t13+t15;
t18 = uh1.^2;
t19 = t6.*t9.*t10.*(1.0./2.0);
t20 = 1.0./t10;
t21 = t12+t13-t15;
t22 = t5.*t16.*5.0e1;
t23 = tanh(t22);
t24 = n1.^2;
t25 = n2.^2;
t26 = t5.*t21.*5.0e1;
t27 = tanh(t26);
t28 = t24.*uh2;
t29 = t25.*uh2;
t30 = n1.*t11.*t14.*uh1.*(1.0./2.0);
t31 = n1.*uh2;
t32 = n2.*uh3;
t33 = t31+t32;
t34 = t5.*t33.*1.0e2;
t35 = tanh(t34);
t36 = param1.*t24.*uh2;
t37 = param1.*t25.*uh2;
t38 = t24.*uh3;
t39 = t25.*uh3;
t40 = n2.*t11.*t14.*uh1.*(1.0./2.0);
t41 = t24+t25;
t42 = param1.*t24.*uh3;
t43 = param1.*t25.*uh3;
t44 = u4-uinf4;
t45 = t4.*t8.*t35.*t41.*uh1.*uh2.*2.0;
t46 = u1-uinf1;
t47 = t28+t29+t30;
t48 = t5.*uh4;
t49 = t5.*t11.*t14.*t16.*t20.*(1.0./4.0);
t50 = t28+t29-t30;
t51 = t5.*t11.*t14.*t20.*t21.*(1.0./4.0);
t52 = -t19+t48+t51;
t53 = param1.*t2;
t54 = param1.*t3;
t69 = param1.*uh1.*uh4.*2.0;
t55 = t2+t3+t53+t54-t69;
t56 = u2-uinf2;
t57 = t28+t29+t30-t36-t37;
t58 = -t28-t29+t30+t36+t37;
t59 = u3-uinf3;
t60 = t41.^2;
t61 = t38+t39+t40-t42-t43;
t62 = -t38-t39+t40+t42+t43;
t63 = t4.*t8.*t35.*t41.*uh1.*uh3.*2.0;
t64 = n1.*uh3;
t72 = n2.*uh2;
t65 = t64-t72;
t66 = t38+t39+t40;
t67 = t19-t48+t49;
t68 = t38+t39-t40;
t70 = n1.*n2.*t35;
t71 = t4.*t8.*t35.*t60.*uh2.*uh3.*2.0;
t73 = t5.*t11.*t14.*t33.*(1.0./2.0);
t74 = t2+t3;
t75 = -t19+t48+t73;
t76 = t19-t48+t73;
t77 = n2.*t5.*t35.*t65;
t78 = n1.*t5.*t35.*t65;
fb = [u1.*(1.0./2.0)-uh1+uinf1.*(1.0./2.0)+t46.*(t4.*t8.*t35.*t55-t4.*t8.*t18.*t27.*t52+t4.*t8.*t18.*t23.*(t19+t49-t5.*uh4)).*(1.0./2.0)-t56.*(t45+t4.*t8.*t20.*t23.*t57.*uh1-t4.*t8.*t20.*t27.*t58.*uh1).*(1.0./2.0)-t59.*(t63+t4.*t8.*t20.*t23.*t61.*uh1-t4.*t8.*t20.*t27.*t62.*uh1).*(1.0./2.0)-t4.*t8.*t18.*t44.*(t23+t27-t35.*2.0).*(1.0./2.0);u2.*(1.0./2.0)-uh2+uinf2.*(1.0./2.0)-t59.*(t70+t71+t4.*t8.*t20.*t23.*t47.*t61-t4.*t8.*t20.*t27.*t50.*t62).*(1.0./2.0)+t56.*(t25.*t35-t2.*t4.*t8.*t35.*t60.*2.0-t4.*t8.*t20.*t23.*t47.*t57+t4.*t8.*t20.*t27.*t50.*t58).*(1.0./2.0)-t44.*(-t45+t4.*t8.*t23.*t47.*uh1+t4.*t8.*t27.*t50.*uh1).*(1.0./2.0)+t46.*(t77-t4.*t8.*t27.*t50.*t52.*uh1+t4.*t8.*t23.*t47.*t67.*uh1+t4.*t5.*t8.*t35.*t41.*t55.*uh2).*(1.0./2.0);u3.*(1.0./2.0)-uh3+uinf3.*(1.0./2.0)-t56.*(t70+t71+t4.*t8.*t20.*t23.*t57.*t66-t4.*t8.*t20.*t27.*t58.*t68).*(1.0./2.0)+t59.*(t24.*t35-t3.*t4.*t8.*t35.*t60.*2.0-t4.*t8.*t20.*t23.*t61.*t66+t4.*t8.*t20.*t27.*t62.*t68).*(1.0./2.0)-t44.*(-t63+t4.*t8.*t23.*t66.*uh1+t4.*t8.*t27.*t68.*uh1).*(1.0./2.0)-t46.*(t78+t4.*t8.*t27.*t52.*t68.*uh1-t4.*t8.*t23.*t66.*t67.*uh1-t4.*t5.*t8.*t35.*t41.*t55.*uh3).*(1.0./2.0);u4.*(1.0./2.0)-uh4+uinf4.*(1.0./2.0)+t46.*(-t9.*t35.*t65.^2+t4.*t8.*t18.*t27.*t52.*t76+t4.*t8.*t18.*t23.*t67.*t75+t4.*t8.*t9.*t35.*t41.*t55.*t74.*(1.0./2.0)).*(1.0./2.0)+t44.*(-t4.*t8.*t18.*t23.*t75+t4.*t8.*t18.*t27.*t76+t4.*t8.*t35.*t41.*t74).*(1.0./2.0)-t56.*(t77+t4.*t5.*t8.*t35.*t60.*t74.*uh2+t4.*t8.*t20.*t23.*t57.*t75.*uh1+t4.*t8.*t20.*t27.*t58.*t76.*uh1).*(1.0./2.0)-t59.*(-t78+t4.*t5.*t8.*t35.*t60.*t74.*uh3+t4.*t8.*t20.*t23.*t61.*t75.*uh1+t4.*t8.*t20.*t27.*t62.*t76.*uh1).*(1.0./2.0)];
