function kkgencodeall(npm, foldername)
% Code generation for all models

kkgencodeelemface("KokkosFlux", npm, 1, foldername);
kkgencodeelemface("KokkosSource", npm, 1, foldername);
kkgencodeelemface("KokkosTdfunc", npm, 1, foldername);

kkgencodeelemface("KokkosEoS", npm, 2, foldername);
kkgencodeelemface("KokkosEoSdu", npm, 2, foldername);
kkgencodeelemface("KokkosEoSdw", npm, 2, foldername);

kkgencodeelemface("KokkosAvfield", npm, 2, foldername);
kkgencodeelemface("KokkosOutput", npm, 2, foldername);
kkgencodeelemface("KokkosSourcew", npm, 2, foldername);

kkgencodeelemface("KokkosFbou", npm, 3, foldername);
kkgencodeelemface("KokkosUbou", npm, 3, foldername);

kkgencodeelemface("KokkosFhat", npm, 4, foldername);
kkgencodeelemface("KokkosUhat", npm, 4, foldername);
kkgencodeelemface("KokkosStab", npm, 4, foldername);

kkgencodeelemface("KokkosInitu", npm, 5, foldername);
kkgencodeelemface("KokkosInitq", npm, 5, foldername);
kkgencodeelemface("KokkosInitudg", npm, 5, foldername);
kkgencodeelemface("KokkosInitwdg", npm, 5, foldername);
kkgencodeelemface("KokkosInitodg", npm, 5, foldername);

kkgencodeelemface("cpuInitu", npm, 5, foldername);
kkgencodeelemface("cpuInitq", npm, 5, foldername);
kkgencodeelemface("cpuInitudg", npm, 5, foldername);
kkgencodeelemface("cpuInitwdg", npm, 5, foldername);
kkgencodeelemface("cpuInitodg", npm, 5, foldername);

kkgencodeelemface("HdgFlux", npm, 6, foldername);
kkgencodeelemface("HdgSource", npm, 6, foldername);
kkgencodeelemface("HdgEoS", npm, 6, foldername);
kkgencodeelemface("HdgSourcew", npm, 6, foldername);
kkgencodeelemface("HdgSourcewonly", npm, 7, foldername);
kkgencodeelemface("HdgFbou", npm, 8, foldername);
kkgencodeelemface("HdgFbouonly", npm, 9, foldername);

end


