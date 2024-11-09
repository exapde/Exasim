function sympyassign2(mystr::String, f, udg, wdg, uhg)

    f = f[:]
    str1 = getccode(f, "f[")
    str1 = "\t\t{\n" * str1 * "\t\t}\n"
    mystr = mystr * str1

    nf = length(f);    
            
    if !isnothing(udg) && !isempty(udg)
        nu = length(udg);
        f_udg = [SymPy.symbols("f_udg$i") for i=1:(nf*nu)];
        for n = 1:nu
          for m = 1:nf      
            f_udg[m + nf*(n-1)] = diff(f[m],udg[n]);      
          end
        end        
        if !isempty(f_udg)
            str2 = getccode(f_udg, "f_udg[")
            str2 = "\t\t{\n" * str2 * "\t\t}\n"
            mystr = mystr * str2
        end
    end

    if !isnothing(wdg) && !isempty(wdg)
        nw = length(wdg);
        f_wdg = [SymPy.symbols("f_wdg$i") for i=1:(nf*nu)];
        for n = 1:nw
          for m = 1:nf      
            f_wdg[m + nf*(n-1)] = diff(f[m],wdg[n]);      
          end
        end                
        if !isempty(f_wdg)
            str3 = getccode(f_wdg, "f_wdg[")
            str3 = "\t\t{\n" * str3 * "\t\t}\n"
            mystr = mystr * str3
        end
    end

    if !isnothing(uhg) && !isempty(uhg)
        nu = length(uhg);
        f_uhg = [SymPy.symbols("f_uhg$i") for i=1:(nf*nu)];
        for n = 1:nu
          for m = 1:nf      
            f_uhg[m + nf*(n-1)] = diff(f[m],uhg[n]);      
          end
        end              
        if !isempty(f_uhg)
            str4 = getccode(f_uhg, "f_uhg[")
            str4 = "\t\t{\n" * str4 * "\t\t}\n"
            mystr = mystr * str4
        end
    end

    return mystr
end

