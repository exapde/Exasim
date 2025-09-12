function gencodeall(npm, foldername)
% Code generation for all models

gencodeelemface("Flux", npm, 1, foldername);
gencodeelemface("Source", npm, 1, foldername);
gencodeelemface("Tdfunc", npm, 1, foldername);

gencodeelemface("EoS", npm, 2, foldername);
gencodeelemface("EoSdu", npm, 2, foldername);
gencodeelemface("EoSdw", npm, 2, foldername);

gencodeelemface("Avfield", npm, 2, foldername);
gencodeelemface("Output", npm, 2, foldername);
gencodeelemface("Sourcew", npm, 2, foldername);

gencodeelemface("Fbou", npm, 3, foldername);
gencodeelemface("Ubou", npm, 3, foldername);

gencodeelemface("Fhat", npm, 4, foldername);
gencodeelemface("Uhat", npm, 4, foldername);
gencodeelemface("Stab", npm, 4, foldername);

gencodeelemface("Initu", npm, 5, foldername);
gencodeelemface("Initq", npm, 5, foldername);
gencodeelemface("Initudg", npm, 5, foldername);
gencodeelemface("Initwdg", npm, 5, foldername);
gencodeelemface("Initodg", npm, 5, foldername);

end

