    program esER
    implicit none
    integer N,i,j,ij,ji,l,Nk,nwmax,mcs,sc,scc,al,kcc
    real(kind=8) pi
    parameter(N=5000,Nk=20,sc=10,pi=4*atan(1.0))
    real(kind=8) opinion(N),self(N),ordersum,coop,cooplevel,alpha,temp
    real(kind=8) beta,nw,x,t,k,cp,kc,p,h,rgre,rgim,rloc(N),rloc_i,cre,cim,rlocn,st,rls,orders,coops,pay(n),paymax
    real(kind=8) opinion_t(N),opinion_rate(N),yk(4,N),A(4),B(4),order,idx(N),orderss,coopss,rlss,z,ri(n),rl,rewardn
    data A/1.0d0,2.0d0,2.0d0,1.0d0/,B/0.5d0,0.50d0,1.0d0,0.0d0/
    integer Nei1(N,NK),Nei2(N,Nk),aaa(N,N),ndeg(N),nei(N,N),Ns(N),S(N),noise(N)
    character(len=100)::filename1,filename2
    h=0.01d0
    beta=1d0
    s=0
    alpha=2.00d0
    write(filename1,'(''alpha='',f4.2)')alpha
    open(10,file="rg\"//trim(adjustl(filename1))//".txt")
    open(11,file="c\"//trim(adjustl(filename1))//".txt")
    open(12,file="rl\"//trim(adjustl(filename1))//".txt")
    do kcc=0,40,1
        kc=kcc*0.025d0
        coopss=0d0
        orderss=0d0
        rlss=0d0
        do scc=1,sc
            call RANDOM_SEED
            call randomnet(nei,ndeg)
            ns=0
            do while(sum(ns)<n/2)
                call random_number(x)
                i=1+floor(x*n)
                Ns(i)=1
            end do
            !!============洗牌
            do i=1,N
                opinion(i)=-pi+i*(2d0*pi)/real(n)
                self(i)=-pi+i*(2d0*pi)/real(n)
            end do
            do i=N,1,-1
                call random_number(x)
                j=1+floor(x*i)
                temp=opinion(i)
                opinion(i)=opinion(j)
                opinion(j)=temp
                call random_number(z)
                j=1+floor(z*i)
                temp=self(i)
                self(i)=self(j)
                self(j)=temp
            end do
            !===============================
            opinion_rate=0d0
            t=0.00d0
            coops=0d0
            orders=0d0
            rls=0d0
            do mcs=1,50000
                t=t+0.01d0
                call RK4(opinion_t,opinion,opinion_rate,h,kc,self,Ns,nei,ndeg)
                rgre=0d0
                rgim=0d0
                pay=0d0
                rl=0d0
                do i=1,N
                    rloc(i)=0d0
                    do ij=1,ndeg(i)
                        j=Nei(i,ij)
                        cre=cos(opinion(i))/2d0+cos(opinion(j))/2d0
                        cim=sin(opinion(i))/2d0+sin(opinion(j))/2d0
                        rloc_i=sqrt((cre)**2d0+(cim)**2d0)
                        rloc(i)=rloc(i)+rloc_i/real(ndeg(i))
                        ri(i)=rloc(i)

                        if (ns(i)==1.and.ns(j)==1)then
                            pay(i)=pay(i)+1d0-((1-rloc_i)**alpha)/2d0
                        elseif (ns(i)==1.and.ns(j)==0)then
                            pay(i)=pay(i)+1d0-(1-rloc_i)**alpha
                        elseif (ns(i)==0.and.ns(j)==1)then
                            pay(i)=pay(i)+1d0
                        elseif (ns(i)==0.and.ns(j)==0)then
                            pay(i)=pay(i)
                        endif

                    end do
                    pay(i)=pay(i)/real(ndeg(i))
                    rgre=rgre+cos(opinion(i))/real(n)
                    rgim=rgim+sin(opinion(i))/real(n)
                    rl=rl+ri(i)/real(n)
                end do
                order=sqrt((rgre)**2d0+(rgim)**2d0)

                !=======================费米规则
                do i=1,N
                    call random_number(x)
                    j=nei(i,1+floor(x*ndeg(i)))
                    call random_number(x)
                    if(x<(1.0d0/(1.0d0+exp(((pay(i)-pay(j)))*beta))))then
                        s(i)=ns(j)
                    else
                        s(i)=ns(i)
                    endif
                enddo
                ns=s

                coop=real(sum(ns))/real(n)
                if (mod(mcs,500)==0)then
                    write(*,'(F5.3, 4X, i2, 4X, i6, 4X, F8.6, 4X, F8.6)')kc,scc,mcs,order,coop
                endif
                if (mcs>40000) then
                    orders=orders+order/10000d0
                    coops=coops+coop/10000d0
                    rls=rls+rl/10000d0
                endif
            end do
            orderss=orderss+orders/real(sc)
            coopss=coopss+coops/real(sc)
            rlss=rlss+rls/real(sc)
        end do
        write(10,"(f12.8)")orderss
        write(11,"(f12.8)")coopss
        write(12,"(f12.8)")rlss
    enddo



    end program

    subroutine fun(opinion_rate,self,opinion,Ns,kc,nei,ndeg)
    integer i,j,NK,N,ij
    parameter(N=5000,NK=20)
    real(kind=8) cp,opinion(N),self(N),opinion_rate(N),kc,rho
    integer nei(N,n),ndeg(N),Ns(N)

    do i=1,N
        cp=0.0d0
        do ij=1,ndeg(i)
            j=Nei(i,ij)
            cp=cp+sin(opinion(j)-opinion(i))
        end do
        opinion_rate(i)=self(i)+kc*cp*Ns(i)
    end do

    return
    end subroutine

    subroutine RK4(opinion_t,opinion,opinion_rate,h,kc,self,Ns,nei,ndeg)
    integer i,j,k,N,NK
    parameter(N=5000,NK=20,pi=4*atan(1.0))
    real(kind=8) p,opinion(N),opinion_rate(N),opinion_t(N),self(N),A(4),B(4),yk(4,N),h,kc,rho
    integer nei(N,n),ndeg(N),Ns(N)
    data A/1.0d0,2.0d0,2.0d0,1.0d0/,B/0.5d0,0.50d0,1.0d0,0.0d0/

    opinion_t=opinion
    do k=1,4
        call fun(opinion_rate,self,opinion_t,Ns,kc,nei,ndeg)
        do i=1,N
            yk(k,i)=opinion_rate(i)
            opinion_t(i)=opinion(i)+h*yk(k,i)*B(k)!y（n+1）
        end do
    end do
    do i=1,N
        p=0.0d0
        do k=1,4
            p=p+(h/6.0d0)*A(k)*yk(k,i)!p指opinion(i)的变化量
        end do
        opinion(i)=opinion(i)+p
        opinion_rate(i)=p/(0.010d0)
        opinion(i) = mod(opinion(i)+pi, 2d0*pi)-pi
    end do

    return
    end subroutine

    subroutine randomnet(nei,ndeg)
    implicit none
    integer N, K_max, K_mean, nedge, nedge2,nedge1, count1,i,j,i1,i2
    parameter (N=5000,K_mean=20,K_max=50,nedge=N*K_mean/2,nedge2=nedge+1000)  !nedge1 为实际边数
    integer nedg(nedge2,2),nei(N,K_max),ndeg(N)
    real(kind=8)z

    !open(11,file="randomnet.dat",status="unknown")
    do i=1,N
        ndeg(i)=0
    enddo

    do i=1,nedge2
        nedg(i,1)=0
        nedg(i,2)=0
    enddo


    do j=1,nedge

10      call random_number(z) !随机找两个合适的点加边
        i1= int(z*N)+1
        if(ndeg(i1).ge.k_max) then
            goto 10
        endif

        count1=0
14      call random_number(z)
        count1=count1+1
        if(count1>nedge) then
            nedge1=j
            goto 12
        endif
        i2= int(z*N)+1
        if(  (ndeg(i2).ge.k_max) .or. (i2.eq.i1)  ) then
            goto 14
        endif

        if(ndeg(i1)>0) then
            do i=1,ndeg(i1)
                if (nei(i1,i).eq.i2) then
                    goto 14
                endif
            enddo
        end if

        nedg(j,1)=i1
        nedg(j,2)=i2
        ndeg(i1)=ndeg(i1)+1
        ndeg(i2)=ndeg(i2)+1

        nei(i1,ndeg(i1))=i2  !添加邻居信息
        nei(i2,ndeg(i2))=i1
        nedge1=j
    enddo

    !==============避免孤立点的存在
12  do i=1,N
        if (ndeg(i).eq.0) then
13          call random_number(z)
            i1= int(z*N)+1
            if(  (ndeg(i1).ge.k_max) .or. (i1.eq.i)  ) then
                goto 13
            end if

            nedge1=nedge1+1
            nedg(nedge1,1)=i
            nedg(nedge1,2)=i1

            ndeg(i)=ndeg(i)+1
            ndeg(i1)=ndeg(i1)+1

            nei(i,ndeg(i))=i1  !添加邻居信息
            nei(i1,ndeg(i1))=i
        endif
    enddo

    return
    end




