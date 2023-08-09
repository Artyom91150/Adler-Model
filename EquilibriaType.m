function N = EquilibriaType(Lambda)

N = 10;

if (imag(Lambda(1)) == 0) || (imag(Lambda(2)) == 0)

    if (abs(real(Lambda(1))) < 1e-1) || (abs(real(Lambda(2))) < 1e-1)
        N = 0; %Сложное СР
        elseif (real(Lambda(1)) * real(Lambda(2)) < 0)
            N = 1; %Седло
            elseif (real(Lambda(1)) < 0) && (real(Lambda(2)) < 0)
                N = 2; %Устойчивый узел
                elseif (real(Lambda(1)) > 0) && (real(Lambda(2)) > 0)
                    N = 3; %Неустойчивый узел
                    
    end
else
    if (abs(real(Lambda(1))) < 1e-5)
        N = 4; %Центр
        elseif (real(Lambda(1)) < 0)
            N = 5; %Устойчивый фокус
            elseif (real(Lambda(1)) > 0)
                N = 6; %Неустойчивый фокус
    end
    
end

end

