    program stability

    implicit none
    integer n,km,k,w,x,i,maxk,fr(201),m,lam,al,ph
    real(kind=8) pi,sigma,mu,a,b,p,g,r,lamuda,r1,r2,r3,e1,alpha,x0,x1,x2,gx,gx1,gx2
    parameter(n=5000,km=20,pi=4*atan(1.0))
    real(kind=8) pk(50),omega(n),e(101),pl,binomial,guass,c1,c2,c3,c4
    character(len=100)::filename0,filename

    pk=0.0d0
    sigma=sqrt(n*(km/real(n-1))*(1-km/real(n-1)))
    mu=km
    pl=km/real(N-1)
    maxk=50
  
    !do k = 0, maxk
    !    pk(k) = binomial(N-1,k)*(pl**k)*((1.0d0-pl)**(N-1-k))
    !    print *, 'Degree k = ', k, ': P(k) = ', pk(k)
    !end do
    do k=1,maxk
        a=k-0.5
        b=k+0.5
        pk(k)=guass(sigma,mu,a,b)
        !open(11,file="pk.txt")
        !write(11,"(f10.8)")pk(k)
    enddo

    print *,sum(pk)
    do w=1,n
        omega(w)=-pi+w*(2d0*pi)/real(n);
    enddo
    
      write(filename,'(''k='',i2,",")') km


    do al=0,100  
        alpha=al*0.02d0
     
   write(filename0,'(''alpha='',f6.3)') alpha
   open(11,file="r\"//trim(adjustl(filename0))//".txt")
   open(12,file="c\"//trim(adjustl(filename0))//".txt")   
    g=0d0
    do x=1,999
         x0=real(x)/1000d0
       gx=sqrt(1d0-x0**2d0)*(1d0-((1d0-sqrt((1d0+x0)/2d0))**alpha)/(2d0-(1d0-sqrt((1d0+x0)/2d0))**alpha))
       if (gx>g) then
           g=gx
       endif
    enddo
    


    do lam=0,100
        lamuda=lam*0.01d0      
        ph=0
        do i=0,100
            r=i*0.01d0
            c4=c3
            r3=0d0
            c3=0d0
            do w=1,n
                r2=0d0
                c2=0d0
                do k=1,maxk
                    if (abs(omega(w))>lamuda*r*k*g) then
                        r1=0d0
                        c1=pk(k)*(1d0-(1d0-2d0/pi)**alpha/(2d0-(1d0-2d0/pi)**alpha))
                    else
                        do x=1,999
                            x0=real(x)/1000d0
                            x1=real(x-1)/1000d0
                            x2=real(x+1)/1000d0
                            !gx1=2d0*sqrt(1d0-x1**2d0)*sqrt((1d0+x1)/2d0)/(1d0+sqrt((1d0+x1)/2d0))
                            !gx2=2d0*sqrt(1d0-x2**2d0)*sqrt((1d0+x2)/2d0)/(1d0+sqrt((1d0+x2)/2d0))
                            gx1=sqrt(1d0-x1**2d0)*(1d0-((1d0-sqrt((1d0+x1)/2d0))**alpha)/(2d0-(1d0-sqrt((1d0+x1)/2d0))**alpha))
                            gx2=sqrt(1d0-x2**2d0)*(1d0-((1d0-sqrt((1d0+x2)/2d0))**alpha)/(2d0-(1d0-sqrt((1d0+x2)/2d0))**alpha))
                            if (gx2<abs(omega(w))/(lamuda*k*r).and.abs(omega(w))/(lamuda*k*r)<gx1) then !ÎÈ¶¨½â
                                exit
                            endif
                        enddo
                        r1=pk(k)*k*x0
                        !c1=pk(k)*2d0*sqrt((1d0+x0)/2d0)/(1d0+sqrt((1d0+x0)/2d0))
                        c1=pk(k)*(1d0-(1d0-sqrt((1d0+x0)/2d0))**alpha/(2d0-(1d0-sqrt((1d0+x0)/2d0))**alpha))
                    end if
                    r2=r2+r1
                    c2=c2+c1
                enddo
                r3=r3+r2/real(n*km)
                c3=c3+c2/real(n)
            enddo
            e(i+1)=r3-r
            !if(e(i+1)<0)then
            !    fr(i+1)=-1
            !elseif (e(i+1)>0)then
            !    fr(i+1)=1
            !else
            !    fr(i+1)=0
            !endif
            write(*,'(F6.4, 2X, F6.4, 2X, F10.7, 2X, F6.4)')lamuda,r,e(i+1),c3
            if (e(i)>0.and.e(i+1)<0)then
            write(11,"(f6.3,4X,f10.8)") lamuda,r-0.005d0
            write(12,"(f6.3,4X,f10.8)") lamuda,(c3+c4)/2d0
            ph=1
            exit
            endif
        enddo      
    enddo
    enddo


    end program stability







    real(kind=8) function binomial(n, k)
    integer  n, k, i


    if (k < 0 .or. k > n) then
        binomial = 0.0d0
    else
        binomial = 1.0d0
        do i = 1, k
            binomial = binomial * (n - i + 1) / i
        end do
    end if
    end function binomial

    real(kind=8) function guass(sigma,mu,a,b) result(p)
    implicit none
    real(kind=8),intent(in)::sigma,mu,a,b
    real(kind=8) pi,h,t
    integer m
    parameter(m=10,pi=4*atan(1.0))
    real(kind=8) xvalue(m+1),yvalue(m+1)
    integer k1,k2,k3

    h=(b-a)/m
    xvalue=0
    yvalue=0
    xvalue(1)=a

    do k1=1,m
        xvalue(k1+1)=xvalue(k1)+h
    enddo

    do k2=1,m+1
        t=xvalue(k2)
        yvalue(k2)=(1/(sigma*sqrt(2*pi)))*exp(-(0.5*(t-mu)**2)/sigma**2)
    enddo

    p=h/3*(yvalue(1)-yvalue(m+1))
    do k3=1,m/2
        p=h/3*(2*yvalue(2*k3+1)+4*yvalue(2*k3))+p
    enddo
    end