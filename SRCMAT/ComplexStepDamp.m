function []=ComplexStepDamp(iPartx,iPartA, x, a)

    [f, fp, fpp]=MathFuncs();
    f2 = @(x,a) f(x+iPartx,a+iPartA);
    fp2 = @(x,a) fp(x+iPartx,a+iPartA);
    fpp2 = @(x,a) fpp(x+iPartx,a+iPartA);
    df = @(x,a) f(x+iPartx,a+iPartA)-f(x,a);
    dfp = @(x,a) fp(x+iPartx,a+iPartA)-fp(x,a);
    dfpp = @(x,a) fpp(x+iPartx,a+iPartA)-fpp(x,a);

    [X,A] = meshgrid(x,a);

    Plots(f2,fp2,fpp2, X, A);
    Plots(df,dfp,dfpp, X, A);
end

function [f, fp, fpp]=MathFuncs()

    f = @(x,a) sqrt(a+x.^2);
    fp = @(x,a) x./sqrt(a+x.^2);
    fpp = @(x,a) 1./(sqrt(a+x.^2).^3);

end

function []=Plots(f,fp,fpp, X, A)
    nP = 3;
    mP = 3;
    h=figure;
    resf =f(X,A);
    resfp=fp(X,A);
    resfpp=fpp(X,A);

    d=@(x) real(x);
    ax(1)=subplot(nP, mP, 1);
    surf(X,A,d(resf));
    ax(2)=subplot(nP, mP, 2);
    surf(X,A,d(resfp));
    ax(3)=subplot(nP, mP, 3);
    surf(X,A,d(resfpp));


    d=@(x) imag(x);
    ax(4)=subplot(nP, mP, 4);
    surf(X,A,d(resf));
    ax(5)=subplot(nP, mP, 5);
    surf(X,A,d(resfp));
    ax(6)=subplot(nP, mP, 6);
    surf(X,A,d(resfpp))


    d=@(x) real(x)+imag(x);
    ax(7)=subplot(nP, mP, 7);
    surf(X,A,d(resf));
    ax(8)=subplot(nP, mP, 8);
    surf(X,A,d(resfp));
    ax(9)=subplot(nP, mP, 9);
    surf(X,A,d(resfpp))

    s=findobj(h,'type','surf');
    [s.LineStyle]  = deal('none');
    [ax.YScale]  = deal('log');
    [ax([3,6]).ZScale]  = deal('log');
end
