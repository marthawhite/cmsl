
if strcmp(computer(), 'PCWIN') == 1 % on windows machine
   
    cmd = 'mex -f C:\Users\xzhang\AppData\Roaming\MathWorks\MATLAB\R2010a\mexopts_ivf.bat -c ';

%    eval([cmd 'tmpbngc.f']);
    eval([cmd 'mpbngc.f']);
    eval([cmd 'plqdf1.f']);
    eval([cmd 'pllpb2.f']);
    eval([cmd 'pbmmex.f']);
    
    mex -output pbm pbm_driver.cpp objfunc.cpp pbmmex.obj mpbngc.obj plqdf1.obj pllpb2.obj
    
else
    
    cmd = 'system (''gfortran -fPIC -c ';
    
    eval([cmd 'pbmmex.f'')']);
    eval([cmd 'mpbngc.f'')']);
    eval([cmd 'plqdf1.f'')']);
    eval([cmd 'pllpb2.f'')']);
    
    mex -lgfortran -output pbm pbm_driver.cpp objfunc.cpp pbmmex.o mpbngc.o plqdf1.o pllpb2.o

end
