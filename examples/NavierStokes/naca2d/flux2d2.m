function [f,f_udg] = flux2d2(pg,udg,param,time)
%FLUX2D2
%    [F,F_UDG] = FLUX2D2(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    29-Mar-2024 17:45:49
[ng,nc] = size(udg);
nch = 4;
nd = 2;
one = ones(ng,1);
param1 = param{1};
param2 = param{2};
param4 = param{4};
param5 = param{5};
param6 = param{6};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
u4 = udg(:,4);
u5 = udg(:,5);
u6 = udg(:,6);
u7 = udg(:,7);
u8 = udg(:,8);
u9 = udg(:,9);
u10 = udg(:,10);
u11 = udg(:,11);
u12 = udg(:,12);
x3 = pg(:,3);
zero = zeros(ng,1);
t2 = u2.^2;
t3 = 1.0./u1.^2;
t4 = t2.*t3.*(1.0./2.0);
t5 = u3.^2;
t6 = t3.*t5.*(1.0./2.0);
t7 = t4+t6;
t14 = t7.*u1;
t8 = -t14+u4;
t9 = param1-1.0;
t10 = 1.0./param1;
t11 = 1.0./param2.^2;
t12 = 1.0./u1;
t13 = 1.0./pi;
t15 = t8.*t9.*3.0e2;
t23 = t10.*t11.*1.5e2;
t16 = t15-t23;
t17 = atan(t16);
t18 = t13.*t17;
t19 = t18+1.0./2.0;
t20 = t8.*t9;
t21 = t10.*t11.*(1.0./2.0);
t22 = param2.^2;
t24 = t20-t21;
t25 = t19.*t24;
t26 = t21+t25+1.061029024220506e-3;
t27 = 1.0./param4;
t28 = param1.*param6.*t12.*t22.*t26;
t29 = t28+5.52e2./5.0;
t30 = 1.0./t29;
t31 = param6+5.52e2./5.0;
t32 = param1.*t12.*t22.*t26;
t33 = t32.^(3.0./2.0);
t41 = t12.*u2.*u5;
t34 = -t41+u6;
t36 = t12.*u3.*u5;
t35 = -t36+u7;
t37 = t12.*t35;
t57 = t12.*u2.*u9;
t38 = -t57+u10;
t39 = t12.*t38;
t40 = t37+t39;
t42 = t12.*t34.*2.0;
t59 = t12.*u3.*u9;
t43 = -t59+u11;
t79 = t12.*t43;
t44 = t42-t79;
t45 = t7.*u5;
t46 = t3.*t34.*u2;
t47 = t3.*t35.*u3;
t48 = t46+t47;
t49 = t48.*u1;
t50 = t45+t49-u8;
t51 = t16.^2;
t52 = t51+1.0;
t53 = 1.0./t52;
t54 = t13.*t16.*t53;
t55 = t18+t54+1.0./2.0;
t56 = t12.*u2.*u3;
t58 = t27.*t30.*t31.*t33.*t40;
t60 = t12.*t26;
t61 = t12.*u4;
t62 = t60+t61;
t63 = t12.*t34;
t113 = t12.*t43.*2.0;
t64 = t63-t113;
t65 = 1.0./param5;
t66 = t7.*u9;
t67 = t3.*t38.*u2;
t68 = t3.*t43.*u3;
t69 = t67+t68;
t70 = t69.*u1;
t71 = t66+t70-u12;
t72 = 1.0./t9;
t73 = 1.0./u1.^3;
t74 = t2.*t73;
t75 = t5.*t73;
t76 = t74+t75;
t78 = t76.*u1;
t77 = t4+t6-t78;
t80 = t9.*t19.*t77;
t81 = t9.*t13.*t24.*t53.*t77.*3.0e2;
t82 = t80+t81;
t83 = param1.*param6.*t3.*t22.*t26;
t84 = param1.*param6.*t12.*t22.*t82;
t85 = t83+t84;
t86 = 1.0./t29.^2;
t87 = param1.*t12.*t22.*t82;
t88 = param1.*t3.*t22.*t26;
t89 = t87+t88;
t90 = sqrt(t32);
t91 = 1.0./u1.^4;
t92 = t3.*t35;
t93 = t3.*t38;
t109 = t73.*u3.*u5;
t110 = t73.*u2.*u9;
t94 = t92+t93-t109-t110;
t95 = t3.*t34.*2.0;
t96 = t73.*u3.*u9;
t162 = t3.*t43;
t97 = t95+t96-t162-t73.*u2.*u5.*2.0;
t98 = t9.*t13.*t53.*t77.*6.0e2;
t99 = 1.0./t52.^2;
t118 = t9.*t13.*t51.*t77.*t99.*6.0e2;
t100 = t98-t118;
t101 = t2.*t91.*u5;
t102 = t5.*t91.*u5;
t103 = t101+t102-t34.*t73.*u2.*2.0-t35.*t73.*u3.*2.0;
t104 = t103.*u1;
t105 = t46+t47+t104-t76.*u5;
t106 = t26.*u5;
t107 = t9.*t50.*t55.*u1;
t108 = t106+t107;
t111 = t27.*t31.*t33.*t40.*t85.*t86;
t112 = t111-t3.*u2.*u3-t27.*t30.*t31.*t33.*t94-t27.*t30.*t31.*t40.*t89.*t90.*(3.0./2.0);
t114 = t3.*t26;
t115 = t3.*u4;
t116 = t12.*t82;
t117 = t114+t115+t116;
t119 = t3.*t34;
t120 = t73.*u3.*u9.*2.0;
t133 = t73.*u2.*u5;
t121 = t119+t120-t133-t3.*t43.*2.0;
t122 = t2.*t91.*u9;
t123 = t5.*t91.*u9;
t124 = t122+t123-t38.*t73.*u2.*2.0-t43.*t73.*u3.*2.0;
t125 = t124.*u1;
t126 = t67+t68+t125-t76.*u9;
t127 = t26.*u9;
t128 = t9.*t55.*t71.*u1;
t129 = t127+t128;
f = [u2+u5.*x3;t21+t2.*t12+u6.*x3+t19.*(t20-t10.*t11.*(1.0./2.0))+t27.*t30.*t31.*t33.*t44.*(2.0./3.0)+1.061029024220506e-3;t56+t58+u7.*x3;t62.*u2+x3.*(u8-t9.*t50.*t55)+t12.*t27.*t30.*t31.*t33.*t40.*u3+t12.*t27.*t30.*t31.*t33.*t44.*u2.*(2.0./3.0)-param1.*t3.*t27.*t30.*t31.*t33.*t65.*t72.*t108;u3+u9.*x3;t56+t58+u10.*x3;t21+t25+t5.*t12+u11.*x3-t27.*t30.*t31.*t33.*t64.*(2.0./3.0)+1.061029024220506e-3;t62.*u3+x3.*(u12-t9.*t55.*t71)+t12.*t27.*t30.*t31.*t33.*t40.*u2-t12.*t27.*t30.*t31.*t33.*t64.*u3.*(2.0./3.0)-param1.*t3.*t27.*t30.*t31.*t33.*t65.*t72.*t129];
if nargout > 1
    t130 = t9.*t12.*t19.*u2;
    t131 = t9.*t12.*t13.*t24.*t53.*u2.*3.0e2;
    t132 = t130+t131;
    t134 = t119-t133;
    t135 = t134.*u1;
    t136 = t3.*u2.*u5;
    t137 = t135+t136;
    t138 = t9.*t12.*t13.*t53.*u2.*6.0e2;
    t144 = t9.*t12.*t13.*t51.*t99.*u2.*6.0e2;
    t139 = t138-t144;
    t140 = param1.^2;
    t141 = t12.*u3;
    t142 = param1.*param6.*t12.*t22.*t27.*t31.*t33.*t40.*t86.*t132;
    t143 = t141+t142-t3.*t27.*t30.*t31.*t33.*u9-param1.*t12.*t22.*t27.*t30.*t31.*t40.*t90.*t132.*(3.0./2.0);
    t145 = t93-t110;
    t146 = t145.*u1;
    t147 = t3.*u2.*u9;
    t148 = t146+t147;
    t149 = t9.*t12.*t19.*u3;
    t150 = t9.*t12.*t13.*t24.*t53.*u3.*3.0e2;
    t151 = t149+t150;
    t152 = t12.*t27.*t30.*t31.*t33.*t40;
    t153 = t92-t109;
    t154 = t153.*u1;
    t155 = t3.*u3.*u5;
    t156 = t154+t155;
    t157 = t9.*t12.*t13.*t53.*u3.*6.0e2;
    t163 = t9.*t12.*t13.*t51.*t99.*u3.*6.0e2;
    t158 = t157-t163;
    t159 = t12.*u2;
    t160 = param1.*param6.*t12.*t22.*t27.*t31.*t33.*t40.*t86.*t151;
    t161 = t159+t160-t3.*t27.*t30.*t31.*t33.*u5-param1.*t12.*t22.*t27.*t30.*t31.*t40.*t90.*t151.*(3.0./2.0);
    t164 = t96-t162;
    t165 = t164.*u1;
    t166 = t9.*t19;
    t167 = param1.*3.0e2;
    t168 = t167-3.0e2;
    t169 = t13.*t24.*t53.*t168;
    t170 = t166+t169;
    t171 = t13.*t53.*t168.*2.0;
    t177 = t13.*t51.*t99.*t168.*2.0;
    t172 = t171-t177;
    t173 = param1.*t12.*t22.*t27.*t30.*t31.*t40.*t90.*t170.*(3.0./2.0);
    t174 = t173-param1.*param6.*t12.*t22.*t27.*t31.*t33.*t40.*t86.*t170;
    t175 = t12.*t170;
    t176 = t12+t175;
    t178 = t12.*t27.*t30.*t31.*t33;
    t179 = t3.*t27.*t30.*t31.*t33.*u2;
    t180 = one.*x3;
    t181 = t9.*t55.*t77.*u1;
    t182 = t21+t25+t181+1.061029024220506e-3;
    t183 = t3.*t27.*t30.*t31.*t33.*u3;
    t184 = t178+x3;
    t185 = t3.*t27.*t30.*t31.*t33.*u2.*(2.0./3.0);
    t186 = t12.*t27.*t30.*t31.*t33.*(4.0./3.0);
    t187 = t186+x3;
    t188 = t9.*t55;
    t189 = t188+1.0;
    t190 = t189.*x3;
    t191 = param1.*t12.*t27.*t30.*t31.*t33.*t55.*t65;
    t192 = t190+t191;
    f_udg = [zero;-t2.*t3-t9.*t19.*t77-t9.*t13.*t24.*t53.*t77.*3.0e2-t27.*t30.*t31.*t33.*t97.*(2.0./3.0)+t27.*t31.*t33.*t44.*t85.*t86.*(2.0./3.0)-t27.*t30.*t31.*t44.*t89.*t90;t112;-t117.*u2+x3.*(t9.*t50.*t100-t9.*t55.*t105)-t3.*t27.*t30.*t31.*t33.*t40.*u3-t3.*t27.*t30.*t31.*t33.*t44.*u2.*(2.0./3.0)-t12.*t27.*t30.*t31.*t33.*t94.*u3-t12.*t27.*t30.*t31.*t33.*t97.*u2.*(2.0./3.0)+t12.*t27.*t31.*t33.*t40.*t85.*t86.*u3+t12.*t27.*t31.*t33.*t44.*t85.*t86.*u2.*(2.0./3.0)-t12.*t27.*t30.*t31.*t40.*t89.*t90.*u3.*(3.0./2.0)-t12.*t27.*t30.*t31.*t44.*t89.*t90.*u2+param1.*t27.*t30.*t31.*t33.*t65.*t72.*t73.*t108.*2.0+param1.*t3.*t27.*t30.*t31.*t33.*t65.*t72.*(t82.*u5-t9.*t50.*t55+t9.*t50.*t100.*u1-t9.*t55.*t105.*u1)-param1.*t3.*t27.*t31.*t33.*t65.*t72.*t85.*t86.*t108+param1.*t3.*t27.*t30.*t31.*t65.*t72.*t89.*t90.*t108.*(3.0./2.0);zero;t112;-t80-t81-t3.*t5+t27.*t30.*t31.*t33.*t121.*(2.0./3.0)-t27.*t31.*t33.*t64.*t85.*t86.*(2.0./3.0)+t27.*t30.*t31.*t64.*t89.*t90;-t117.*u3+x3.*(t9.*t71.*t100-t9.*t55.*t126)-t3.*t27.*t30.*t31.*t33.*t40.*u2+t3.*t27.*t30.*t31.*t33.*t64.*u3.*(2.0./3.0)-t12.*t27.*t30.*t31.*t33.*t94.*u2+t12.*t27.*t30.*t31.*t33.*t121.*u3.*(2.0./3.0)+t12.*t27.*t31.*t33.*t40.*t85.*t86.*u2-t12.*t27.*t30.*t31.*t40.*t89.*t90.*u2.*(3.0./2.0)-t12.*t27.*t31.*t33.*t64.*t85.*t86.*u3.*(2.0./3.0)+t12.*t27.*t30.*t31.*t64.*t89.*t90.*u3+param1.*t27.*t30.*t31.*t33.*t65.*t72.*t73.*t129.*2.0+param1.*t3.*t27.*t30.*t31.*t33.*t65.*t72.*(t82.*u9-t9.*t55.*t71+t9.*t71.*t100.*u1-t9.*t55.*t126.*u1)-param1.*t3.*t27.*t31.*t33.*t65.*t72.*t85.*t86.*t129+param1.*t3.*t27.*t30.*t31.*t65.*t72.*t89.*t90.*t129.*(3.0./2.0);one;t12.*u2.*2.0-t9.*t12.*t19.*u2-t9.*t12.*t13.*t24.*t53.*u2.*3.0e2-t3.*t27.*t30.*t31.*t33.*u5.*(4.0./3.0)-param1.*t12.*t22.*t27.*t30.*t31.*t44.*t90.*t132+param1.*param6.*t12.*t22.*t27.*t31.*t33.*t44.*t86.*t132.*(2.0./3.0);t143;t60+t61+x3.*(t9.*t50.*t139-t9.*t55.*t137)-t12.*t132.*u2+t12.*t27.*t30.*t31.*t33.*t44.*(2.0./3.0)-t27.*t30.*t31.*t33.*t73.*u2.*u5.*(4.0./3.0)-t27.*t30.*t31.*t33.*t73.*u3.*u9+param1.*t3.*t27.*t30.*t31.*t33.*t65.*t72.*(t132.*u5+t9.*t50.*t139.*u1-t9.*t55.*t137.*u1)-param1.*t3.*t22.*t27.*t30.*t31.*t40.*t90.*t132.*u3.*(3.0./2.0)-param1.*t3.*t22.*t27.*t30.*t31.*t44.*t90.*t132.*u2+param1.*param6.*t3.*t22.*t27.*t31.*t33.*t40.*t86.*t132.*u3+param1.*param6.*t3.*t22.*t27.*t31.*t33.*t44.*t86.*t132.*u2.*(2.0./3.0)+t22.*t27.*t30.*t31.*t65.*t72.*t73.*t90.*t108.*t132.*t140.*(3.0./2.0)-param6.*t22.*t27.*t31.*t33.*t65.*t72.*t73.*t86.*t108.*t132.*t140;zero;t143;-t130-t131+t3.*t27.*t30.*t31.*t33.*u5.*(2.0./3.0)+param1.*t12.*t22.*t27.*t30.*t31.*t64.*t90.*t132-param1.*param6.*t12.*t22.*t27.*t31.*t33.*t64.*t86.*t132.*(2.0./3.0);t152-x3.*(t9.*t55.*t148-t9.*t71.*t139)-t12.*t132.*u3+t27.*t30.*t31.*t33.*t73.*u3.*u5.*(2.0./3.0)-t27.*t30.*t31.*t33.*t73.*u2.*u9+param1.*t3.*t27.*t30.*t31.*t33.*t65.*t72.*(t132.*u9-t9.*t55.*t148.*u1+t9.*t71.*t139.*u1)-param1.*t3.*t22.*t27.*t30.*t31.*t40.*t90.*t132.*u2.*(3.0./2.0)+param1.*t3.*t22.*t27.*t30.*t31.*t64.*t90.*t132.*u3+param1.*param6.*t3.*t22.*t27.*t31.*t33.*t40.*t86.*t132.*u2-param1.*param6.*t3.*t22.*t27.*t31.*t33.*t64.*t86.*t132.*u3.*(2.0./3.0)+t22.*t27.*t30.*t31.*t65.*t72.*t73.*t90.*t129.*t132.*t140.*(3.0./2.0)-param6.*t22.*t27.*t31.*t33.*t65.*t72.*t73.*t86.*t129.*t132.*t140;zero;-t9.*t12.*t19.*u3-t9.*t12.*t13.*t24.*t53.*u3.*3.0e2+t3.*t27.*t30.*t31.*t33.*u9.*(2.0./3.0)-param1.*t12.*t22.*t27.*t30.*t31.*t44.*t90.*t151+param1.*param6.*t12.*t22.*t27.*t31.*t33.*t44.*t86.*t151.*(2.0./3.0);t161;t152+x3.*(t9.*t50.*t158-t9.*t55.*t156)-t12.*t151.*u2-t27.*t30.*t31.*t33.*t73.*u3.*u5+t27.*t30.*t31.*t33.*t73.*u2.*u9.*(2.0./3.0)+param1.*t3.*t27.*t30.*t31.*t33.*t65.*t72.*(t151.*u5+t9.*t50.*t158.*u1-t9.*t55.*t156.*u1)-param1.*t3.*t22.*t27.*t30.*t31.*t40.*t90.*t151.*u3.*(3.0./2.0)-param1.*t3.*t22.*t27.*t30.*t31.*t44.*t90.*t151.*u2+param1.*param6.*t3.*t22.*t27.*t31.*t33.*t40.*t86.*t151.*u3+param1.*param6.*t3.*t22.*t27.*t31.*t33.*t44.*t86.*t151.*u2.*(2.0./3.0)+t22.*t27.*t30.*t31.*t65.*t72.*t73.*t90.*t108.*t140.*t151.*(3.0./2.0)-param6.*t22.*t27.*t31.*t33.*t65.*t72.*t73.*t86.*t108.*t140.*t151;one;t161;-t149-t150+t12.*u3.*2.0-t3.*t27.*t30.*t31.*t33.*u9.*(4.0./3.0)+param1.*t12.*t22.*t27.*t30.*t31.*t64.*t90.*t151-param1.*param6.*t12.*t22.*t27.*t31.*t33.*t64.*t86.*t151.*(2.0./3.0);t60+t61+x3.*(t9.*t71.*t158+t9.*t55.*(t165-t3.*u3.*u9))-t12.*t151.*u3-t12.*t27.*t30.*t31.*t33.*t64.*(2.0./3.0)-t27.*t30.*t31.*t33.*t73.*u2.*u5-t27.*t30.*t31.*t33.*t73.*u3.*u9.*(4.0./3.0)+param1.*t3.*t27.*t30.*t31.*t33.*t65.*t72.*(t151.*u9+t9.*t71.*t158.*u1+t9.*t55.*u1.*(t165-t3.*u3.*u9))-param1.*t3.*t22.*t27.*t30.*t31.*t40.*t90.*t151.*u2.*(3.0./2.0)+param1.*t3.*t22.*t27.*t30.*t31.*t64.*t90.*t151.*u3+param1.*param6.*t3.*t22.*t27.*t31.*t33.*t40.*t86.*t151.*u2-param1.*param6.*t3.*t22.*t27.*t31.*t33.*t64.*t86.*t151.*u3.*(2.0./3.0)+t22.*t27.*t30.*t31.*t65.*t72.*t73.*t90.*t129.*t140.*t151.*(3.0./2.0)-param6.*t22.*t27.*t31.*t33.*t65.*t72.*t73.*t86.*t129.*t140.*t151;zero;t166+t169+param1.*t12.*t22.*t27.*t30.*t31.*t44.*t90.*t170-param1.*param6.*t12.*t22.*t27.*t31.*t33.*t44.*t86.*t170.*(2.0./3.0);t174;t176.*u2-t9.*t50.*t172.*x3-param1.*t3.*t27.*t30.*t31.*t33.*t65.*t72.*(t170.*u5+t9.*t50.*t172.*u1)+param1.*t3.*t22.*t27.*t30.*t31.*t40.*t90.*t170.*u3.*(3.0./2.0)+param1.*t3.*t22.*t27.*t30.*t31.*t44.*t90.*t170.*u2-param1.*param6.*t3.*t22.*t27.*t31.*t33.*t40.*t86.*t170.*u3-param1.*param6.*t3.*t22.*t27.*t31.*t33.*t44.*t86.*t170.*u2.*(2.0./3.0)-t22.*t27.*t30.*t31.*t65.*t72.*t73.*t90.*t108.*t140.*t170.*(3.0./2.0)+param6.*t22.*t27.*t31.*t33.*t65.*t72.*t73.*t86.*t108.*t140.*t170;zero;t174;t166+t169-param1.*t12.*t22.*t27.*t30.*t31.*t64.*t90.*t170+param1.*param6.*t12.*t22.*t27.*t31.*t33.*t64.*t86.*t170.*(2.0./3.0);t176.*u3-t9.*t71.*t172.*x3-param1.*t3.*t27.*t30.*t31.*t33.*t65.*t72.*(t170.*u9+t9.*t71.*t172.*u1)+param1.*t3.*t22.*t27.*t30.*t31.*t40.*t90.*t170.*u2.*(3.0./2.0)-param1.*t3.*t22.*t27.*t30.*t31.*t64.*t90.*t170.*u3-param1.*param6.*t3.*t22.*t27.*t31.*t33.*t40.*t86.*t170.*u2+param1.*param6.*t3.*t22.*t27.*t31.*t33.*t64.*t86.*t170.*u3.*(2.0./3.0)-t22.*t27.*t30.*t31.*t65.*t72.*t73.*t90.*t129.*t140.*t170.*(3.0./2.0)+param6.*t22.*t27.*t31.*t33.*t65.*t72.*t73.*t86.*t129.*t140.*t170;t180;t3.*t27.*t30.*t31.*t33.*u2.*(-4.0./3.0);-t3.*t27.*t30.*t31.*t33.*u3;-t9.*t55.*t77.*x3-t2.*t27.*t30.*t31.*t33.*t73.*(4.0./3.0)-t5.*t27.*t30.*t31.*t33.*t73-param1.*t3.*t27.*t30.*t31.*t33.*t65.*t72.*t182;zero;-t3.*t27.*t30.*t31.*t33.*u3;t185;t27.*t30.*t31.*t33.*t73.*u2.*u3.*(-1.0./3.0);zero;t187;zero;-t9.*t12.*t55.*u2.*x3+t3.*t27.*t30.*t31.*t33.*u2.*(4.0./3.0)-param1.*t3.*t27.*t30.*t31.*t33.*t55.*t65.*u2;zero;zero;t12.*t27.*t30.*t31.*t33.*(-2.0./3.0);t3.*t27.*t30.*t31.*t33.*u3.*(-2.0./3.0);zero;zero;t184;t183-t9.*t12.*t55.*u3.*x3-param1.*t3.*t27.*t30.*t31.*t33.*t55.*t65.*u3;zero;t178;zero;t179;zero;zero;zero;t192;zero;zero;zero;zero;zero;t3.*t27.*t30.*t31.*t33.*u3.*(2.0./3.0);-t179;t27.*t30.*t31.*t33.*t73.*u2.*u3.*(-1.0./3.0);t180;-t179;t3.*t27.*t30.*t31.*t33.*u3.*(-4.0./3.0);-t9.*t55.*t77.*x3-t2.*t27.*t30.*t31.*t33.*t73-t5.*t27.*t30.*t31.*t33.*t73.*(4.0./3.0)-param1.*t3.*t27.*t30.*t31.*t33.*t65.*t72.*t182;zero;zero;t178;t183;zero;t184;zero;t179-t9.*t12.*t55.*u2.*x3-param1.*t3.*t27.*t30.*t31.*t33.*t55.*t65.*u2;zero;t12.*t27.*t30.*t31.*t33.*(-2.0./3.0);zero;-t185;zero;zero;t187;-t9.*t12.*t55.*u3.*x3+t3.*t27.*t30.*t31.*t33.*u3.*(4.0./3.0)-param1.*t3.*t27.*t30.*t31.*t33.*t55.*t65.*u3;zero;zero;zero;zero;zero;zero;zero;t192];
end
f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);