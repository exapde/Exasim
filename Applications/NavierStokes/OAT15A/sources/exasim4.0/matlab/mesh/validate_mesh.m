function validate_mesh(dmdmat,dmdcpp)

Compare(dmdmat.nbsd, 1+dmdcpp.nbsd, 'nbsd: ');
Compare(dmdmat.intelem, 1+dmdcpp.intelem, 'intelem: ');
Compare(dmdmat.intent, 1+dmdcpp.intent, 'intent: ');
Compare(dmdmat.elempart, 1+dmdcpp.elempart, 'elempart: ');
Compare(dmdmat.elempartpts, dmdcpp.elempartpts, 'elempartpts: ');
Compare(dmdmat.entpart, 1+dmdcpp.entpart, 'entpart: ');
Compare(dmdmat.entpartpts, dmdcpp.entpartpts, 'entpartpts: ');
Compare(dmdmat.elemrecv, 1+dmdcpp.elemrecv, 'elemrecv: ');
Compare(dmdmat.elemrecvpts, dmdcpp.elemrecvpts, 'elemrecvpts: ');
Compare(dmdmat.elemsend, 1+dmdcpp.elemsend, 'elemsend: ');
Compare(dmdmat.elemsendpts, dmdcpp.elemsendpts, 'elemsendpts: ');        
Compare(dmdmat.entrecv, 1+dmdcpp.entrecv, 'entrecv: ');
Compare(dmdmat.entrecvpts, dmdcpp.entrecvpts, 'entrecvpts: ');
Compare(dmdmat.entsend, 1+dmdcpp.entsend, 'entsend: ');
Compare(dmdmat.entsendpts, dmdcpp.entsendpts, 'entsendpts: ');
Compare(dmdmat.vecrecv, 1+dmdcpp.vecrecv, 'vecrecv: ');
Compare(dmdmat.vecrecvpts, dmdcpp.vecrecvpts, 'vecrecvpts: ');
Compare(dmdmat.vecsend, 1+dmdcpp.vecsend, 'vecsend: ');
Compare(dmdmat.vecsendpts, dmdcpp.vecsendpts, 'vecsendpts: ');
Compare(dmdmat.matrecv, dmdcpp.matrecv, 'matrecv: ');
Compare(dmdmat.matrecvpts, dmdcpp.matrecvpts, 'matrecvpts: ');
Compare(dmdmat.matsend, dmdcpp.matsend, 'matsend: ');
Compare(dmdmat.matsendpts, dmdcpp.matsendpts, 'matsendpts: ');
Compare(dmdmat.rowent2elem, dmdcpp.rowent2elem, 'rowent2elem: ');
Compare(dmdmat.colent2elem, 1+dmdcpp.colent2elem, 'colent2elem: ');
Compare(dmdmat.rowent2ent, dmdcpp.rowent2ent, 'rowent2ent: ');
Compare(dmdmat.colent2ent, 1+dmdcpp.colent2ent, 'colent2ent: ');
Compare(dmdmat.bcrs_rowent2elem, dmdcpp.bcrs_rowent2elem, 'bcrs_rowent2elem: ');
Compare(dmdmat.bcrs_colent2elem, 1+dmdcpp.bcrs_colent2elem, 'bcrs_colent2elem: ');
Compare(dmdmat.bcrs_rowent2ent, dmdcpp.bcrs_rowent2ent, 'bcrs_rowent2ent: ');
Compare(dmdmat.bcrs_colent2ent, 1+dmdcpp.bcrs_colent2ent, 'bcrs_colent2ent: ');        
Compare(dmdmat.ent2ind, 1+dmdcpp.ent2ind, 'ent2ind: ');
Compare(dmdmat.elcon, 1+dmdcpp.elcon, 'elcon: ');
Compare(dmdmat.t2f, 1+dmdcpp.t2f, 't2f: ');

function Compare(A, B, str)

if numel(A) ~= numel(B)    
    disp([str 'A and B have a different number of elements']);
    disp(['size A: ' num2str(size(A))]);
    disp(['size B: ' num2str(size(B))]);
else
    e = max(abs(A(:)-B(:)));
    disp([str 'Maximum absolute error between A and B is ' num2str(e)]);    
end












