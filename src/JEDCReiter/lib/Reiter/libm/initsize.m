function z = initsize(varargin)


    isderiv = 0;
    for i=1:length(varargin)
        if (isa(varargin{i},'deriv1'))
            isderiv=1;
            %== number of variable to take derivative ==%
            nder = nindep(varargin{i});

            %== break from for? ==%
            %what if I have more than one?
            break;
        end
    end

    %== z is of size varargin{1} ==%
    if (isderiv==0)
        z = 0*varargin{1};  % automatically creates sparse matrix if_ input is sparse;
    else
        z = deriv1(zeros(size(varargin{1})),0,nder);
    end
