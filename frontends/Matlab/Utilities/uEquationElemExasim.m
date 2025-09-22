function [Ru, B, D, fg, fg_udg, sg, sg_udg, xg, udg] = uEquationElemExasim(dataout, pde, nge, npe, ne)

%npe = size(mesh.xpe,1);
%ne = size(mesh.t, 2);
nc = pde.nc;
ncu = pde.ncu;
ncq = pde.ncq;
nd = pde.nd;

filename = dataout + "uEquationElem_B.bin";
tmp = readbin(filename);
B = reshape(tmp, npe, npe, ne, ncu, ncq);

filename = dataout + "uEquationElem_D.bin";
tmp = readbin(filename);
D = reshape(tmp, npe, npe, ne, ncu, ncu);

filename = dataout + "uEquationElem_Ru.bin";
tmp = readbin(filename);
Ru = reshape(tmp, npe, ne, ncu);

if nargout>3
  filename = dataout + "uEquationElem_sg.bin";
  tmp = readbin(filename);
  sg = reshape(tmp, nge*ne, ncu);

  filename = dataout + "uEquationElem_sg_udg.bin";
  tmp = readbin(filename);
  sg_udg = reshape(tmp, nge*ne, ncu, nc);

  filename = dataout + "uEquationElem_fg.bin";
  tmp = readbin(filename);
  fg = reshape(tmp, nge*ne, ncu, nd);

  filename = dataout + "uEquationElem_fg_udg.bin";
  tmp = readbin(filename);
  fg_udg = reshape(tmp, nge*ne, ncu, nd, nc);

  filename = dataout + "uEquationElem_xg.bin";
  tmp = readbin(filename);
  xg = reshape(tmp, nge*ne, nd);

  filename = dataout + "uEquationElem_udg.bin";
  tmp = readbin(filename);
  udg = reshape(tmp, nge*ne, nc);
else
  fg = [];
  fg_udg = [];
  sg = [];
  sg_udg = [];
  xg = [];
  udg = [];
end

