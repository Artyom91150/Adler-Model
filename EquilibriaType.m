function N = EquilibriaType(Lambda)

N = 10;

if (imag(Lambda(1)) == 0) || (imag(Lambda(2)) == 0)

    if (abs(real(Lambda(1))) < 1e-1) || (abs(real(Lambda(2))) < 1e-1)
        N = 0; %������� ��
        elseif (real(Lambda(1)) * real(Lambda(2)) < 0)
            N = 1; %�����
            elseif (real(Lambda(1)) < 0) && (real(Lambda(2)) < 0)
                N = 2; %���������� ����
                elseif (real(Lambda(1)) > 0) && (real(Lambda(2)) > 0)
                    N = 3; %������������ ����
                    
    end
else
    if (abs(real(Lambda(1))) < 1e-5)
        N = 4; %�����
        elseif (real(Lambda(1)) < 0)
            N = 5; %���������� �����
            elseif (real(Lambda(1)) > 0)
                N = 6; %������������ �����
    end
    
end

end

