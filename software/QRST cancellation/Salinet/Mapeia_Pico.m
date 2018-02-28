function[Pico] = Mapeia_Pico(ref,sinal)

global Fam

TamJanela = round(0.010*Fam);
if abs(sinal(ref)) > abs(sinal(ref-1)) & abs(sinal(ref)) > abs(sinal(ref+1))
    Pico = ref;
else
    
    index = ref;
    ref1 = ref;
    while index >= max(1,ref1-TamJanela)
        if abs(sinal(index)) > abs(sinal(ref))
            ref = index;
        end
        index = index - 1;
    end


    index = ref1;
    ref1 = ref;
    ref2 = index;
    ref = index;
    
    while index <= min(length(sinal),ref2+TamJanela)
        if abs(sinal(index)) > abs(sinal(ref))
            ref = index;
        end
        index = index + 1;
    end
    
    ref2 = ref;
    if abs(sinal(ref1)) > abs(sinal(ref2))
        Pico = ref1;
    else Pico = ref2;
    end

end


    
    
            
    
    