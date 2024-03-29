function gencodeall(npm)
% Code generation for all models

gencodeelemface("Flux", npm, 1);
gencodeelemface("Source", npm, 1);
gencodeelemface("Tdfunc", npm, 1);

gencodeelemface("EoS", npm, 2);
gencodeelemface("EoSdu", npm, 2);
gencodeelemface("EoSdw", npm, 2);

gencodeelemface("Avfield", npm, 2);
gencodeelemface("Output", npm, 2);
gencodeelemface("Sourcew", npm, 2);

gencodeelemface("Fbou", npm, 3);
gencodeelemface("Ubou", npm, 3);

gencodeelemface("Fhat", npm, 4);
gencodeelemface("Uhat", npm, 4);
gencodeelemface("Stab", npm, 4);

gencodeelemface("Initu", npm, 5);
gencodeelemface("Initq", npm, 5);
gencodeelemface("Initudg", npm, 5);
gencodeelemface("Initwdg", npm, 5);
gencodeelemface("Initodg", npm, 5);

end

