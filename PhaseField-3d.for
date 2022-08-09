
      module kvisual
      implicit none
      real*8 UserVar(70000,13,10)
      integer nelem
      save
      end module
      
      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     3 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njpro,period)

      use kvisual
      include 'aba_param.inc' !implicit real(a-h o-z)
      
      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*),
     1 energy(*),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel),
     2 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)

      parameter(ndim=3,ntens=6,nsvint=14)
      
      dimension dN(nnode,1),dNdz(ndim,nnode),stran(ntens),
     2 dNdx(ndim,nnode),b(ntens,nnode*ndim),ddsdde(ntens,ntens),
     3 stress(ntens),dstran(ntens),statevLocal(nsvint)
      
      
!     initialising
      do k1=1,ndofel
       rhs(k1,1)=0.d0
      end do
      amatrx=0.d0
      
!     find number of elements          
      if (dtime.eq.0.d0) then
       if (jelem.eq.1) then
        nelem=jelem
       else
        if (jelem.gt.nelem) nelem=jelem 
       endif 
      endif      
      
!     reading parameters
      xlc=props(3)
      Gc=props(4)
      xk=props(5)
      
!     C3D4      
      if(jtype .eq .1)then
      wght=0.166666667
      ninpt=4
      do kintk=1,ninpt
!     evaluate shape functions and derivatives
       call kshapefcn1(kintk,ninpt,nnode,ndim,dN,dNdz)      
       call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
       dvol=wght*djac
       
!     form B-matrix
       b=0.d0
       do inod=1,nnode
        b(1,3*inod-2)=dNdx(1,inod)
        b(2,3*inod-1)=dNdx(2,inod)
        b(3,3*inod)=  dNdx(3,inod)
        b(4,3*inod-2)=dNdx(2,inod)
        b(4,3*inod-1)=dNdx(1,inod)
        b(5,3*inod-1)=dNdx(3,inod)
        b(5,3*inod)=  dNdx(2,inod)
        b(6,3*inod-2)=dNdx(3,inod)
        b(6,3*inod)=  dNdx(1,inod)
       end do                     
       
!     compute from nodal values
       phi=0.d0
       do inod=1,nnode
        phi=phi+dN(inod,1)*u(ndim*nnode+inod)
       end do   
       if (phi.gt.1.d0) phi=1.d0
           
!     compute the increment of strain and recover history variables
       dstran=matmul(b,du(1:ndim*nnode,1))
       
       call kstatevar(kintk,nsvint,svars,statevLocal,1) 
       
       stress=statevLocal(1:ntens)
       stran(1:ntens)=statevLocal((ntens+1):(2*ntens))
       phin=statevLocal(2*ntens+1)
       Hn=statevLocal(2*ntens+2)
       if (dtime.eq.0.d0) phin=phi
       
!     call umat to obtain stresses and constitutive matrix 
       call kumat(props,ddsdde,stress,dstran,ntens,statevLocal)
       stran=stran+dstran       
       
!     compute strain energy density from the previous increment       
       Psi=0.d0
       do k1=1,ntens
        Psi=Psi+stress(k1)*stran(k1)*0.5d0
       end do     
       
!     enforcing Karush-Kuhn-Tucker conditions
       if (Psi.gt.Hn) then
        H=Psi
       else
        H=Hn
       endif
       
       statevLocal(1:ntens)=stress(1:ntens)
       statevLocal((ntens+1):(2*ntens))=stran(1:ntens)
       statevLocal(2*ntens+1)=phi
       statevLocal(2*ntens+2)=H
      
       call kstatevar(kintk,nsvint,svars,statevLocal,0)
       
       amatrx(1:12,1:12)=amatrx(1:12,1:12)+
     1 dvol*(((1.d0-phi)**2+xk)*matmul(matmul(transpose(b),ddsdde),b))
        
       rhs(1:12,1)=rhs(1:12,1)-
     1 dvol*(matmul(transpose(b),stress)*((1.d0-phi)**2+xk))       
           
        amatrx(13:16,13:16)=amatrx(13:16,13:16)
     1 +dvol*(matmul(transpose(dNdx),dNdx)*Gc*xlc
     2 +matmul(dN,transpose(dN))*(Gc/xlc+2.d0*H))   

        rhs(13:16,1)=rhs(13:16,1)
     1 -dvol*(matmul(transpose(dNdx),matmul(dNdx,u(13:16)))
     2 *Gc*xlc+dN(:,1)*((Gc/xlc+2.d0*H)*phi-2.d0*H))           
                   
! output
       UserVar(jelem,1:6,kintk)=statevLocal(1:ntens)*((1.d0-phi)**2+xk)
       UserVar(jelem,7:13,kintk)=statevLocal((ntens+1):(2*ntens+1))
      
      end do       ! end loop on material integration points
! C3D6
      elseif(jtype .eq .2)then
      wght=0.166666667
      ninpt=6
      do kintk=1,ninpt
!     evaluate shape functions and derivatives
       call kshapefcn2(kintk,ninpt,nnode,ndim,dN,dNdz)      
       call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
       dvol=wght*djac
       
!     form B-matrix
       b=0.d0
       do inod=1,nnode
        b(1,3*inod-2)=dNdx(1,inod)
        b(2,3*inod-1)=dNdx(2,inod)
        b(3,3*inod)=  dNdx(3,inod)
        b(4,3*inod-2)=dNdx(2,inod)
        b(4,3*inod-1)=dNdx(1,inod)
        b(5,3*inod-1)=dNdx(3,inod)
        b(5,3*inod)=  dNdx(2,inod)
        b(6,3*inod-2)=dNdx(3,inod)
        b(6,3*inod)=  dNdx(1,inod)
       end do                     
       
!     compute from nodal values
       phi=0.d0
       do inod=1,nnode
        phi=phi+dN(inod,1)*u(ndim*nnode+inod)
       end do   
       if (phi.gt.1.d0) phi=1.d0
           
!     compute the increment of strain and recover history variables
       dstran=matmul(b,du(1:ndim*nnode,1))
       
       call kstatevar(kintk,nsvint,svars,statevLocal,1) 
       
       stress=statevLocal(1:ntens)
       stran(1:ntens)=statevLocal((ntens+1):(2*ntens))
       phin=statevLocal(2*ntens+1)
       Hn=statevLocal(2*ntens+2)
       if (dtime.eq.0.d0) phin=phi
       
!     call umat to obtain stresses and constitutive matrix 
       call kumat(props,ddsdde,stress,dstran,ntens,statevLocal)
       stran=stran+dstran
       
!     compute strain energy density from the previous increment       
       Psi=0.d0
       do k1=1,ntens
        Psi=Psi+stress(k1)*stran(k1)*0.5d0
       end do     

!     enforcing Karush-Kuhn-Tucker conditions
       if (Psi.gt.Hn) then
        H=Psi
       else
        H=Hn
       endif
       
       statevLocal(1:ntens)=stress(1:ntens)
       statevLocal((ntens+1):(2*ntens))=stran(1:ntens)
       statevLocal(2*ntens+1)=phi
       statevLocal(2*ntens+2)=H
      
       call kstatevar(kintk,nsvint,svars,statevLocal,0)
       
       amatrx(1:18,1:18)=amatrx(1:18,1:18)+
     1 dvol*(((1.d0-phi)**2+xk)*matmul(matmul(transpose(b),ddsdde),b))
        
       rhs(1:18,1)=rhs(1:18,1)-
     1 dvol*(matmul(transpose(b),stress)*((1.d0-phi)**2+xk))       
           
        amatrx(19:24,19:24)=amatrx(19:24,19:24)
     1 +dvol*(matmul(transpose(dNdx),dNdx)*Gc*xlc
     2 +matmul(dN,transpose(dN))*(Gc/xlc+2.d0*H))   

        rhs(19:24,1)=rhs(19:24,1)
     1 -dvol*(matmul(transpose(dNdx),matmul(dNdx,u(19:24)))
     2 *Gc*xlc+dN(:,1)*((Gc/xlc+2.d0*H)*phi-2.d0*H))           
                   
! output
       UserVar(jelem,1:6,kintk)=statevLocal(1:ntens)*((1.d0-phi)**2+xk)
       UserVar(jelem,7:13,kintk)=statevLocal((ntens+1):(2*ntens+1))
      
      end do       ! end loop on material integration points
! C3D8
      elseif(jtype .eq .3)then
      wght=1.d0
      ninpt=8
      do kintk=1,ninpt
!     evaluate shape functions and derivatives
       call kshapefcn3(kintk,ninpt,nnode,ndim,dN,dNdz)      
       call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
       dvol=wght*djac
       
!     form B-matrix
       b=0.d0
       do inod=1,nnode
        b(1,3*inod-2)=dNdx(1,inod)
        b(2,3*inod-1)=dNdx(2,inod)
        b(3,3*inod)=  dNdx(3,inod)
        b(4,3*inod-2)=dNdx(2,inod)
        b(4,3*inod-1)=dNdx(1,inod)
        b(5,3*inod-1)=dNdx(3,inod)
        b(5,3*inod)=  dNdx(2,inod)
        b(6,3*inod-2)=dNdx(3,inod)
        b(6,3*inod)=  dNdx(1,inod)
       end do                     
       
!     compute from nodal values
       phi=0.d0
       do inod=1,nnode
        phi=phi+dN(inod,1)*u(ndim*nnode+inod)
       end do   
       if (phi.gt.1.d0) phi=1.d0
           
!     compute the increment of strain and recover history variables
       dstran=matmul(b,du(1:ndim*nnode,1))
       
       call kstatevar(kintk,nsvint,svars,statevLocal,1) 
       
       stress=statevLocal(1:ntens)
       stran(1:ntens)=statevLocal((ntens+1):(2*ntens))
       phin=statevLocal(2*ntens+1)
       Hn=statevLocal(2*ntens+2)
       if (dtime.eq.0.d0) phin=phi
       
!     call umat to obtain stresses and constitutive matrix 
       call kumat(props,ddsdde,stress,dstran,ntens,statevLocal)
       stran=stran+dstran
       
!     compute strain energy density from the previous increment       
       Psi=0.d0
       do k1=1,ntens
        Psi=Psi+stress(k1)*stran(k1)*0.5d0
       end do     

!     enforcing Karush-Kuhn-Tucker conditions
       if (Psi.gt.Hn) then
        H=Psi
       else
        H=Hn
       endif
       
       statevLocal(1:ntens)=stress(1:ntens)
       statevLocal((ntens+1):(2*ntens))=stran(1:ntens)
       statevLocal(2*ntens+1)=phi
       statevLocal(2*ntens+2)=H
      
       call kstatevar(kintk,nsvint,svars,statevLocal,0)
       
       amatrx(1:24,1:24)=amatrx(1:24,1:24)+
     1 dvol*(((1.d0-phi)**2+xk)*matmul(matmul(transpose(b),ddsdde),b))
        
       rhs(1:24,1)=rhs(1:24,1)-
     1 dvol*(matmul(transpose(b),stress)*((1.d0-phi)**2+xk))       
           
        amatrx(25:32,25:32)=amatrx(25:32,25:32)
     1 +dvol*(matmul(transpose(dNdx),dNdx)*Gc*xlc
     2 +matmul(dN,transpose(dN))*(Gc/xlc+2.d0*H))   

        rhs(25:32,1)=rhs(25:32,1)
     1 -dvol*(matmul(transpose(dNdx),matmul(dNdx,u(25:32)))
     2 *Gc*xlc+dN(:,1)*((Gc/xlc+2.d0*H)*phi-2.d0*H))           
                   
! output
       UserVar(jelem,1:6,kintk)=statevLocal(1:ntens)*((1.d0-phi)**2+xk)
       UserVar(jelem,7:13,kintk)=statevLocal((ntens+1):(2*ntens+1))
      
      end do       ! end loop on material integration points
! C3D10
      elseif(jtype .eq. 4)then
      wght=0.041666667
      ninpt=10
      do kintk=1,ninpt
!     evaluate shape functions and derivatives
       call kshapefcn4(kintk,ninpt,nnode,ndim,dN,dNdz)      
       call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
       dvol=wght*djac
       
!     form B-matrix
       b=0.d0
       do inod=1,nnode
        b(1,3*inod-2)=dNdx(1,inod)
        b(2,3*inod-1)=dNdx(2,inod)
        b(3,3*inod)=  dNdx(3,inod)
        b(4,3*inod-2)=dNdx(2,inod)
        b(4,3*inod-1)=dNdx(1,inod)
        b(5,3*inod-1)=dNdx(3,inod)
        b(5,3*inod)=  dNdx(2,inod)
        b(6,3*inod-2)=dNdx(3,inod)
        b(6,3*inod)=  dNdx(1,inod)
       end do                     

!     compute from nodal values
       phi=0.d0
       do inod=1,nnode
        phi=phi+dN(inod,1)*u(ndim*nnode+inod)
       end do   
       if (phi.gt.1.d0) phi=1.d0

!     compute the increment of strain and recover history variables
       dstran=matmul(b,du(1:ndim*nnode,1))
       
       call kstatevar(kintk,nsvint,svars,statevLocal,1) 
       
       stress=statevLocal(1:ntens)
       stran(1:ntens)=statevLocal((ntens+1):(2*ntens))
       phin=statevLocal(2*ntens+1)
       Hn=statevLocal(2*ntens+2)
       if (dtime.eq.0.d0) phin=phi
       
!     call umat to obtain stresses and constitutive matrix 
       call kumat(props,ddsdde,stress,dstran,ntens,statevLocal)
       stran=stran+dstran
       
!     compute strain energy density from the previous increment       
       Psi=0.d0
       do k1=1,ntens
        Psi=Psi+stress(k1)*stran(k1)*0.5d0
       end do     

!     enforcing Karush-Kuhn-Tucker conditions
       if (Psi.gt.Hn) then
        H=Psi
       else
        H=Hn
       endif
       
       statevLocal(1:ntens)=stress(1:ntens)
       statevLocal((ntens+1):(2*ntens))=stran(1:ntens)
       statevLocal(2*ntens+1)=phi
       statevLocal(2*ntens+2)=H
      
       call kstatevar(kintk,nsvint,svars,statevLocal,0)
       
       amatrx(1:30,1:30)=amatrx(1:30,1:30)+
     1 dvol*(((1.d0-phi)**2+xk)*matmul(matmul(transpose(b),ddsdde),b))
        
       rhs(1:30,1)=rhs(1:30,1)-
     1 dvol*(matmul(transpose(b),stress)*((1.d0-phi)**2+xk))       
           
        amatrx(31:40,31:40)=amatrx(31:40,31:40)
     1 +dvol*(matmul(transpose(dNdx),dNdx)*Gc*xlc
     2 +matmul(dN,transpose(dN))*(Gc/xlc+2.d0*H))   

        rhs(31:40,1)=rhs(31:40,1)
     1 -dvol*(matmul(transpose(dNdx),matmul(dNdx,u(31:40)))
     2 *Gc*xlc+dN(:,1)*((Gc/xlc+2.d0*H)*phi-2.d0*H))           
                   
! output
       UserVar(jelem,1:6,kintk)=statevLocal(1:ntens)*((1.d0-phi)**2+xk)
       UserVar(jelem,7:13,kintk)=statevLocal((ntens+1):(2*ntens+1))
      
      end do       ! end loop on material integration points
      endif
      RETURN
      END
c*****************************************************************      
      subroutine kshapefcn1(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      dimension dN(nnode,1),dNdz(ndim,*),coord34(3,4)
      
      data  coord34 /0.13819660,0.13819660,0.13819660,
     2               0.58541020,0.13819660,0.13819660,
     3               0.13819660,0.58541020,0.13819660,
     4               0.13819660,0.13819660,0.58541020/       
     
!     3D 4-nodes

!     determine (g,h,r)
      g=coord34(1,kintk)
      h=coord34(2,kintk)
      r=coord34(3,kintk)

!     shape functions 
      dN(1,1)=1.d0-g-h-r
      dN(2,1)=g
      dN(3,1)=h
      dN(4,1)=r

!     derivative d(Ni)/d(g)
      dNdz(1,1)=-1.d0
      dNdz(1,2)=1.d0
      dNdz(1,3)=0.d0
      dNdz(1,4)=0.d0

!     derivative d(Ni)/d(h)
      dNdz(2,1)=-1.d0
      dNdz(2,2)=0.d0
      dNdz(2,3)=1.d0
      dNdz(2,4)=0.d0
      
!     derivative d(Ni)/d(r)
      dNdz(3,1)=-1.d0
      dNdz(3,2)=0.d0
      dNdz(3,3)=0.d0
      dNdz(3,4)=1.d0
      
      return
      end 
c*****************************************************************
      subroutine kshapefcn2(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      dimension dN(nnode,1),dNdz(ndim,*),coord36(3,6)
      
      data  coord36 / 0.166666667,0.166666667,-0.577350269,
     2                0.666666667,0.166666667,-0.577350269,
     3                0.166666667,0.666666667,-0.577350269,
     4                0.166666667,0.166666667, 0.577350269,
     5                0.666666667,0.166666667, 0.577350269,
     6                0.166666667,0.666666667, 0.577350269/      
     
!     3D 6-nodes

!     determine (g,h,r)
      g=coord36(1,kintk)
      h=coord36(2,kintk)
      r=coord36(3,kintk)

!     shape functions 
      dN(1,1)=(1.d0-g-h)*(1.d0-r)/2.d0
      dN(2,1)=g*(1.d0-r)/2.d0
      dN(3,1)=h*(1.d0-r)/2.d0
      dN(4,1)=(1.d0-g-h)*(1.d0+r)/2.d0
      dN(5,1)=g*(1.d0+r)/2.d0
      dN(6,1)=h*(1.d0+r)/2.d0   
   
      
!     derivative d(Ni)/d(g)
      dNdz(1,1)=-(1.d0-r)/2.d0
      dNdz(1,2)= (1.d0-r)/2.d0
      dNdz(1,3)=0.d0
      dNdz(1,4)=-(1.d0+r)/2.d0
      dNdz(1,5)= (1.d0+r)/2.d0
      dNdz(1,6)=0.d0

!     derivative d(Ni)/d(h)
      dNdz(2,1)=-(1.d0-r)/2.d0
      dNdz(2,2)=0.d0
      dNdz(2,3)= (1.d0-r)/2.d0
      dNdz(2,4)=-(1.d0+r)/2.d0
      dNdz(2,5)=0.d0
      dNdz(2,6)= (1.d0+r)/2.d0
      
!     derivative d(Ni)/d(r)
      dNdz(3,1)=-(1.d0-g-h)/2.d0
      dNdz(3,2)=-g/2.d0
      dNdz(3,3)=-h/2.d0
      dNdz(3,4)= (1.d0-g-h)/2.d0
      dNdz(3,5)= g/2.d0
      dNdz(3,6)= h/2.d0
      
      return
      end 
c*****************************************************************
      subroutine kshapefcn3(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      dimension dN(nnode,1),dNdz(ndim,*),coord38(3,8)
      
      data  coord38 /-0.577350269,-0.577350269,-0.577350269,
     2                0.577350269,-0.577350269,-0.577350269,
     3               -0.577350269, 0.577350269,-0.577350269,
     4                0.577350269, 0.577350269,-0.577350269,
     5               -0.577350269,-0.577350269, 0.577350269,
     6                0.577350269,-0.577350269, 0.577350269,
     7               -0.577350269, 0.577350269, 0.577350269,
     8                0.577350269, 0.577350269, 0.577350269/      
     
!     3D 8-nodes

!     determine (g,h,r)
      g=coord38(1,kintk)
      h=coord38(2,kintk)
      r=coord38(3,kintk)

!     shape functions 
      dN(1,1)=((1.d0-g)*(1.d0-h)*(1.d0-r))/8.d0
      dN(2,1)=((1.d0+g)*(1.d0-h)*(1.d0-r))/8.d0
      dN(3,1)=((1.d0+g)*(1.d0+h)*(1.d0-r))/8.d0
      dN(4,1)=((1.d0-g)*(1.d0+h)*(1.d0-r))/8.d0
      dN(5,1)=((1.d0-g)*(1.d0-h)*(1.d0+r))/8.d0
      dN(6,1)=((1.d0+g)*(1.d0-h)*(1.d0+r))/8.d0     
      dN(7,1)=((1.d0+g)*(1.d0+h)*(1.d0+r))/8.d0        
      dN(8,1)=((1.d0-g)*(1.d0+h)*(1.d0+r))/8.d0        
      
!     derivative d(Ni)/d(g)
      dNdz(1,1)=-(1.d0-h)*(1.d0-r)/8.d0
      dNdz(1,2)= (1.d0-h)*(1.d0-r)/8.d0
      dNdz(1,3)= (1.d0+h)*(1.d0-r)/8.d0
      dNdz(1,4)=-(1.d0+h)*(1.d0-r)/8.d0
      dNdz(1,5)=-(1.d0-h)*(1.d0+r)/8.d0
      dNdz(1,6)= (1.d0-h)*(1.d0+r)/8.d0
      dNdz(1,7)= (1.d0+h)*(1.d0+r)/8.d0 
      dNdz(1,8)=-(1.d0+h)*(1.d0+r)/8.d0   

!     derivative d(Ni)/d(h)
      dNdz(2,1)=-(1.d0-g)*(1.d0-r)/8.d0
      dNdz(2,2)=-(1.d0+g)*(1.d0-r)/8.d0
      dNdz(2,3)= (1.d0+g)*(1.d0-r)/8.d0
      dNdz(2,4)= (1.d0-g)*(1.d0-r)/8.d0
      dNdz(2,5)=-(1.d0-g)*(1.d0+r)/8.d0
      dNdz(2,6)=-(1.d0+g)*(1.d0+r)/8.d0
      dNdz(2,7)= (1.d0+g)*(1.d0+r)/8.d0
      dNdz(2,8)= (1.d0-g)*(1.d0+r)/8.d0
      
!     derivative d(Ni)/d(r)
      dNdz(3,1)=-(1.d0-g)*(1.d0-h)/8.d0
      dNdz(3,2)=-(1.d0+g)*(1.d0-h)/8.d0
      dNdz(3,3)=-(1.d0+g)*(1.d0+h)/8.d0
      dNdz(3,4)=-(1.d0-g)*(1.d0+h)/8.d0
      dNdz(3,5)= (1.d0-g)*(1.d0-h)/8.d0
      dNdz(3,6)= (1.d0+g)*(1.d0-h)/8.d0
      dNdz(3,7)= (1.d0+g)*(1.d0+h)/8.d0
      dNdz(3,8)= (1.d0-g)*(1.d0+h)/8.d0
      
      return
      end 
c*****************************************************************
      subroutine kshapefcn4(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      dimension dN(nnode,1),dNdz(ndim,*),coord34(3,4)
      
      data  coord34 /0.13819660,0.13819660,0.13819660,
     2               0.58541020,0.13819660,0.13819660,
     3               0.13819660,0.58541020,0.13819660,
     4               0.13819660,0.13819660,0.58541020/     
     
!     3D 10-nodes

!     determine (g,h,r)
      g=coord34(1,kintk)
      h=coord34(2,kintk)
      r=coord34(3,kintk)
      
!     shape functions 
      dN(1,1)=(2.d0*(1.d0-g-h-r)-1.d0)*(1.d0-g-h-r)
      dN(2,1)=(2.d0*g-1.d0)*g
      dN(3,1)=(2.d0*h-1.d0)*h
      dN(4,1)=(2.d0*r-1.d0)*r
      dN(5,1)=4.d0*(1.d0-g-h-r)*g
      dN(6,1)=4.d0*g*h    
      dN(7,1)=4.d0*(1.d0-g-h-r)*h       
      dN(8,1)=4.d0*(1.d0-g-h-r)*r
      dN(9,1)=4.d0*g*r     
      dN(10,1)=4.d0*h*r   
    
!     derivative d(Ni)/d(g)
      dNdz(1,1)=4.d0*(g+h+r)-3.d0
      dNdz(1,2)=4.d0*g-1.d0
      dNdz(1,3)=0.d0
      dNdz(1,4)=0.d0
      dNdz(1,5)=4.d0*(1-2.d0*g-h-r)
      dNdz(1,6)=4.d0*h
      dNdz(1,7)=-4.d0*h 
      dNdz(1,8)=-4.d0*r
      dNdz(1,9)=4.d0*r  
      dNdz(1,10)=0.d0  

!     derivative d(Ni)/d(h)
      dNdz(2,1)=4.d0*(g+h+r)-3.d0
      dNdz(2,2)=0.d0
      dNdz(2,3)=4.d0*h-1.d0
      dNdz(2,4)=0.d0
      dNdz(2,5)=-4.d0*g
      dNdz(2,6)=4.d0*g
      dNdz(2,7)=4.d0*(1-g-2.d0*h-r)
      dNdz(2,8)=-4.d0*r
      dNdz(2,9)=0.d0      
      dNdz(2,10)=4.d0*r
      
!     derivative d(Ni)/d(r)
      dNdz(3,1)=4.d0*(g+h+r)-3.d0
      dNdz(3,2)=0.d0
      dNdz(3,3)=0.d0
      dNdz(3,4)=4.d0*r-1.d0
      dNdz(3,5)=-4.d0*g
      dNdz(3,6)=0.d0
      dNdz(3,7)=-4.d0*h
      dNdz(3,8)=4.d0*(1-g-h-2.d0*r)
      dNdz(3,9)=4.d0*g     
      dNdz(3,10)=4.d0*h
      
      return
      end 
c*****************************************************************
      subroutine kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
!     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix 
!     dNdx - shape functions derivatives w.r.t. global coordinates
      include 'aba_param.inc'

      dimension xjac(ndim,ndim),xjaci(ndim,ndim),coords(mcrd,nnode),
     1 dNdz(ndim,nnode),dNdx(ndim,nnode)

      xjac=0.d0

      do inod=1,nnode
       do idim=1,ndim
        do jdim=1,ndim
         xjac(jdim,idim)=xjac(jdim,idim)+
     1        dNdz(jdim,inod)*coords(idim,inod)      
        end do
       end do 
      end do

      djac=xjac(1,1)*xjac(2,2)*xjac(3,3)+xjac(1,2)*xjac(2,3)*xjac(3,1)+
     1  xjac(2,1)*xjac(3,2)*xjac(1,3)-xjac(1,3)*xjac(2,2)*xjac(3,1)-
     2  xjac(1,2)*xjac(2,1)*xjac(3,3)-xjac(2,3)*xjac(3,2)*xjac(1,1)
      if (djac.gt.0.d0) then ! jacobian is positive - o.k.
       xjaci(1,1)= (xjac(2,2)*xjac(3,3)-xjac(2,3)*xjac(3,2))/djac
       xjaci(1,2)=-(xjac(1,2)*xjac(3,3)-xjac(3,2)*xjac(1,3))/djac
       xjaci(1,3)= (xjac(1,2)*xjac(2,3)-xjac(1,3)*xjac(2,2))/djac
       xjaci(2,1)=-(xjac(2,1)*xjac(3,3)-xjac(2,3)*xjac(3,1))/djac
       xjaci(2,2)= (xjac(1,1)*xjac(3,3)-xjac(1,3)*xjac(3,1))/djac       
       xjaci(2,3)=-(xjac(1,1)*xjac(2,3)-xjac(1,3)*xjac(2,1))/djac   
       xjaci(3,1)= (xjac(2,1)*xjac(3,2)-xjac(2,2)*xjac(3,1))/djac     
       xjaci(3,2)=-(xjac(1,1)*xjac(3,2)-xjac(1,2)*xjac(3,1))/djac       
       xjaci(3,3)= (xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1))/djac 
         
      else ! negative or zero jacobian
       write(*,*)'WARNING: element',jelem,'has neg. Jacobian'
      endif
	  
      dNdx=matmul(xjaci,dNdz) 
			
      return
      end

c*****************************************************************
      subroutine kstatevar(npt,nsvint,statev,statev_ip,icopy)
c
c     Transfer data to/from element-level state variable array from/to
c     material-point level state variable array.
c
      include 'aba_param.inc'

      dimension statev(*),statev_ip(*)

      isvinc=(npt-1)*nsvint     ! integration point increment

      if (icopy.eq.1) then ! Prepare arrays for entry into umat
       do i=1,nsvint
        statev_ip(i)=statev(i+isvinc)
       enddo
      else ! Update element state variables upon return from umat
       do i=1,nsvint
        statev(i+isvinc)=statev_ip(i)
       enddo
      end if

      return
      end

c*****************************************************************
      subroutine kumat(props,ddsdde,stress,dstran,ntens,statev)
c
c     Subroutine with the material model
c
      include 'aba_param.inc' !implicit real(a-h o-z)
      
      dimension props(*),ddsdde(ntens,ntens),stress(ntens),statev(*),
     + dstran(ntens)

!     Initialization
      ddsdde=0.d0
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      
!     Build stiffness matrix
      eg2=E/(1.d0+xnu)
      elam=(E/(1.d0-2.d0*xnu)-eg2)/3.d0
      
!     Update stresses
      do k1=1,3
       do k2=1,3
        ddsdde(k2,k1)=elam
       end do
       ddsdde(k1,k1)=eg2+elam
      end do
      ddsdde(4,4)=eg2/2.d0
      ddsdde(5,5)=eg2/2.d0
      ddsdde(6,6)=eg2/2.d0
      
      stress=stress+matmul(ddsdde,dstran)   

      return
      end

c*****************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use kvisual
      include 'aba_param.inc' !implicit real(a-h o-z)

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3),jstep(4)

      ddsdde=0.0d0
      noffset=noel-nelem
      statev(1:nstatv)=UserVar(noffset,1:nstatv,npt)
     
      return
      end
