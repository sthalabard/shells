!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            SHELL ROUTINES IN FORTRAN  TO SPEED UP THE RK INTEGRATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            POUR COMPILER : f2py -c -m Shell4KH_boost Shell4KH_boost.f90 --fcompiler=gnu95 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!LES TABLEAUX DE COMPLEXES SONT REPRESENTES PAR DES REELS
! POUR DES RAISONS MYSTERIEUSES, F2PY SEMBLE BUGGER SINON 

!! RHS NON-LINEAIRE
SUBROUTINE Compute_DE(DE,u,v,n)
    implicit none
    integer,intent(in):: n
    real(8),intent(in),dimension(0:n-1,0:1):: u,v
    real(8),intent(out):: DE
    integer::k
    
    DE=0.
    do k=0,(n-1)
        DE= DE+0.5*(v(k,0)-u(k,0))**2+(v(k,1)-u(k,1))**2
    end do
END SUBROUTINE

!! RHS NON-LINEAIRE
SUBROUTINE RHS(da,a,gam,lam,Etas,n)
    implicit none
    integer,intent(in):: n
    real(8),intent(in):: gam,lam
    real(8),intent(in),dimension(0:n-1):: Etas
    real(8),intent(in),dimension(0:n-1,0:1):: a
    real(8),intent(out),dimension(0:n-1,0:1):: da
    real(8):: coeff
    integer:: k,km,kp

        do k =0,(n-1)
                da(k,0)=0.
                da(k,1)=0.
                if (k<n-2) then  !a(k+1)*a(k+2)
                        km=k+1
                        kp=k+2
                        coeff=(lam**2)*gam
                        da(k,0) =da(k,0)+ coeff * (a(km,0)*a(kp,0)-a(km,1)*a(kp,1))
                        da(k,1) =da(k,1)+ coeff * (-a(km,0)*a(kp,1)-a(km,1)*a(kp,0))
                endif

                if ((k<n-1) .AND. (k>0)) then !a(k-1)*a(k+1)
                        km=k-1
                        kp=k+1
                        coeff=(-lam)*(1.+gam)
                        da(k,0) =da(k,0)+ coeff * (a(km,0)*a(kp,0)-a(km,1)*a(kp,1))
                        da(k,1) =da(k,1)+ coeff * (-a(km,0)*a(kp,1)-a(km,1)*a(kp,0))
                endif

                if (k >1) then !a(k-2)a(k-1)
                        km=k-2
                        kp=k-1
                        coeff=1.
                        da(k,0) =da(k,0)+ coeff * (a(km,0)*a(kp,0)-a(km,1)*a(kp,1))
                        da(k,1) =da(k,1)+ coeff * (-a(km,0)*a(kp,1)-a(km,1)*a(kp,0))
                endif

                da(k,0)=da(k,0)*Etas(k) 
                da(k,1)=da(k,1)*Etas(k)
        end do
END SUBROUTINE
 
 
!! UPDATE VIA RUNGE KUTTA
SUBROUTINE UPDATE(u,v,dt,nt,tol,dissip,gam,lam,Etas,n,nstep)
        implicit none
        integer,intent(in):: n,nt
        integer,intent(out):: nstep
        !f2py intent(out) :: nstep
        real(8),intent(in):: dt,gam,lam,tol
        real(8),intent(in),dimension(0:n-1):: dissip,Etas
        real(8),intent(inout),dimension(0:n-1,0:1):: u,v
        !f2py intent(in,out,inplace) :: u,v

        real(8),dimension(0:n-1,0:1):: du,dv,u1,u2,u3,u4,v1,v2,v3,v4
        integer:: i
        real(8)::DE
        
        nstep=0
        do i=1,nt
            call Compute_DE(DE,u,v,n)
            if (DE>tol) then
                 EXIT
            endif
            nstep=nstep+1
            !transfer nl u 
            call rhs(u1,u,gam,lam,etas,n)
            call rhs(u2,u+0.5*dt*u1,gam,lam,etas,n)
            call rhs(u3,u+0.5*dt*u2,gam,lam,etas,n)
            call rhs(u4,u+ dt*u3,gam,lam,etas,n)
            du=(u1+2.*u2+ 2.*u3+u4)/6.
            !dissip
            u(:,0)=(u(:,0)+dt*du(:,0))*exp(-dissip(:)*dt)
            u(:,1)=(u(:,1)+dt*du(:,1))*exp(-dissip(:)*dt)

            !transfer nl v 
            call rhs(v1,v,gam,lam,etas,n)
            call rhs(v2,v+0.5*dt*v1,gam,lam,etas,n)
            call rhs(v3,v+0.5*dt*v2,gam,lam,etas,n)
            call rhs(v4,v+ dt*v3,gam,lam,etas,n)
            dv=(v1+2.*v2+ 2.*v3+v4)/6.
            !dissip
            v(:,0)=(v(:,0)+dt*dv(:,0))*exp(-dissip(:)*dt)
            v(:,1)=(v(:,1)+dt*dv(:,1))*exp(-dissip(:)*dt)

        end do

END SUBROUTINE


