 
 !#include 'Abaqus_Definitions.for'      
 !DIR$ FREEFORM
 module kvisual
 use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 implicit none 
 REAL(dp):: USERVAR(157734,7,4), Ac 
 REAL(dp):: E_11, E_22, MU_12, NU_12, G_F,G_MI, G_MII, ALPHA, BETA
 REAL(dp):: THETA, NVEC(2),LB,OMEGA_MAT(2,2)
 REAL(dp):: STIFF(3,3),Q_MAT(3,3)
 REAL(dp):: w_qp(4),qpT(2,4),qp(4,2)
 real(dp):: ndstor(157734,4,4) !nstor(nelem,ninpt,nnode)
 real(dp):: det_jacbstor(157734,4,1)!jacstor(nelem,ninpt,1)
 real(dp):: bdstor(157734,4,2,4) !bdstor(nelem,ninpt,ndim,nnode)
 real(dp):: Tbar(3,3), Tbar_trans(3,3)
 real(dp):: alphT 
 integer :: nelemp!,nelemc
 save
 end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars, &
 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime, &
 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf, &
 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njpro,period)

 use kvisual
 use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind

 implicit none 
 integer :: nrhs, nsvars, nprops, mcrd, nnode, jtype, kstep, kinc, jelem
 integer :: ndload, npredf, mlvarx, mdload, njpro, ndofel 
 
 real(dp) :: rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*), & 
 energy(*),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel), &
 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*), &
 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*) 
 real(dp) :: period, pnewdt, dtime 

 if(maxval(u) .gt. 1.0d9) then 
     write(6,*) 'the max of u for element',jelem, 'tends to inf', &
    'hence, terminating the analysis here!'
      call xit 
 endif 
 call element_isotropic(rhs,amatrx,props,coords,u,svars, &
 jelem,dtime,kinc, nsvars, du, kstep)  

 end subroutine UEL   
!================subroutine element_isotropic========================================================     
!!!           in: props,coords,u,svars
!!!           out: rhs,amatrx
 subroutine element_isotropic(rhs,amatrx,props,coords,u,svars,&
            jelem,dtime,kinc, nsvars, du, kstep)
 use kvisual
  use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 implicit none
!coordinate 1、2、3、4
!          x1 x2 x3 x4
!          y1 y2 y3 y4      
 real(dp):: rhs(12,1), amatrx(12,12),props(10), coords(2,4), u(12)
 real(dp):: svars(24), dtime, du(12,1), ddd(4,1)
 integer:: jelem,nsvars,kinc 
!local varibles
 real(dp):: ku(8,8),kd(4,4),ru(8,1),rd(4,1),uu(8,1),dd(4,1)      
 integer:: i, j, k, K1, k2, k3
 real(dp):: det_jacb
 real(dp):: n(4), nd(1,4), b(3,8), bd(2,4) 
 real(dp):: ppf(1,1), pf, dppf(1,1), dpf 
 real(dp):: dalpha, ddalpha, c0      
 real(dp):: omega, domega, ddomega 
 real(dp):: ph_energy, d_ph_energy            
 real(dp):: kd1(4,4), rd1(4,1)
 real(dp):: nd_t(4,1), bd_t(4,2), b_t(8,3)
 real(dp):: bd_t_mul_bd(4,4), bd_t_mul_bd_mul_dd(4,1)
 real(dp):: nd_t_mul_nd(4,4), b_t_mul_d(8,3)
 real(dp):: b_t_mul_d_mul_b(8,8)
 real(dp):: psip
 real(dp):: xi,eta
 integer:: incvar 
 real(dp):: b_t_mul_d_mul_b_mul_uu(8,1)
 real(dp):: dd_copy(4,1)
 REAL(dp):: BD_T_MUL_OMEGA(4,2)
 real(dp):: alphn, alph, alphBn, alphB, Fdeg!, Hfn, Hf, HmIn, HmI, HmIIn, HmII 
 ! real(dp):: pfn, omegan, domegan, ddomegan 
 ! real(dp):: psi_f, psi_mI, psi_mII
 real(dp):: coordx 
 REAL(dp), parameter:: PI = 3.1415926535897932384626433832_dp
 real(dp):: delta_alphn, delta_alph 
 real(dp):: pfn, omegan, domegan, ddomegan, psipn  
 integer:: kstep 
 ! real(dp):: ddse, ddsen !damage degraded strain energy 
 real(dp):: act_psip, act_psipn !actual positive part of strain energy 
 real(dp):: ratio !ratio of omega/omegan 
 real(dp), parameter:: cycles = 1.0_dp 
! if(kinc .gt. 10) call xit 

!INITIALIZE THE KVISUAL VARIABLES ONLY FOR THE FIRST INCREMENT
IF(DTIME .EQ. 0.0_dp) THEN 
   DO I=1,24
      SVARS(I)=0.0_dp
   ENDDO 

!   find number of elements (phase-field)         
   IF(JELEM .EQ. 1) THEN 
       NELEMP = JELEM 
   ELSE 
     IF(JELEM .GT. NELEMP) NELEMP = JELEM
   ENDIF    

!   INITIALIZATION OF MATERIAL PARAMTERS 
   E_11 = PROPS(1) !LONGITUDINAL MODULUS 
   E_22 = PROPS(2) !TRANSVERSE MODULUS
   MU_12 = PROPS(3) !INPLANE SHEAR-MODULUS 
   NU_12 = PROPS(4) !INPLANE POISSON'S RATIO
   G_F = PROPS(5) !FIBER FRACTURE TOUGHNESS 
   G_MI = PROPS(6) !MATRIX MODE-I FRACTURE TOUGHNESS 
   G_MII = PROPS(7) !MATRIX MODE-II FRACTURE TOUGHNESS 
   LB = PROPS(8) !LENGTH SCALE PARAMETER 
   BETA = PROPS(9) !MODEL PENALTY PARAMETER 
   ALPHA = PROPS(10) !ANGLE(IN DEGREES) OF FIBER-ORIENTATION (WRT GLOBAL X-Y)

   alphT = 2.30_dp !Fatigue threshold 
   
!   CALCULATE THE ANGLE ALPHA IN RADIANS 
   THETA = (ALPHA*PI)/180.0_dp  
   
!   CALCULATE THE N VECTOR 
   NVEC(1) = DCOS(THETA)
   NVEC(2) = DSIN(THETA)
!    NVEC(1) = 1.0_dp/SQRT(2.0_dp); NVEC(2) = 1.0_dp/SQRT(2.0_dp)
   
!   CALCULATE THE OMEGA MATRIX 
   DO K1=1,2 
      DO K2=1,2 
         IF(K1 .EQ. K2) THEN 
             OMEGA_MAT(K1,K2) = 1.0_dp + BETA*NVEC(K1)*NVEC(K2)
         ELSE 
             OMEGA_MAT(K1,K2) = BETA*NVEC(K1)*NVEC(K2)
         ENDIF
      ENDDO 
   ENDDO
   
!  CALCULATE THE STIFFNESS AND COMPLIANCE MATRICES 
   CALL CALCULATE_STIFF_COMP(STIFF,PROPS,JELEM)
   
!  CALCULATE THE TRANSFORMATION MATRIX 
   CALL CALCULATE_TRANS(THETA,JELEM,Tbar,Tbar_trans)
   
!  CALCULATE THE CONSTITUTIVE MATRIX 
   CALL Q_MATRIX(Q_MAT,Tbar_trans,STIFF,Tbar,JELEM) 
      
!  Evaluate the gauss weights & co-ordinates 
   CALL GET_W_QP(W_QP,QPT,QP) 

   do k=1,4
      !Get xi and eta 
      xi = qp(k,1)
      eta = qp(k,2) 
      call fem(xi,eta,coords,jelem,k) 
   enddo 
ENDIF 

! if(kstep .eq. 2) then
!    amatrx = 0.0_dp  
!    do i =1,12 
!       amatrx(i,i) = 1.0_dp 
!       rhs(i,1) = 0.0_dp
!    enddo 
! endif 

!initialize varibles
 ku=0.0_dp; kd=0.0_dp; rd=0.0_dp; ru=0.0_dp 
 
 do i=1,4
     uu(2*i-1,1)=u(3*i-2)
     uu(2*i,1)=u(3*i-1)
     dd(i,1)=u(3*i)
     ddd(i,1) = du(3*i,1)
     dd_copy(i,1)=u(3*i)
 enddo

 do k = 1, 4

   nd(1,1:4) = ndstor(jelem,k,1:4)
   det_jacb = det_jacbstor(jelem,k,1)
   bd = 0.0_dp
   bd(1,1:4) = bdstor(jelem,k,1,1:4)
   bd(2,1:4) = bdstor(jelem,k,2,1:4)
   b = 0.0_dp 
   b(1,1:7:2) = bd(1,1:4)
   b(2,2:8:2) = bd(2,1:4)
   b(3,1:7:2) = bd(2,1:4)
   b(3,2:8:2) = bd(1,1:4)

   incvar = (k-1)*6
   psipn = svars(incvar+1)
   pfn = svars(incvar+2)
   alphn = svars(incvar+3)
   alphBn = svars(incvar+4)
   delta_alphn = svars(incvar+5)
   act_psipn = svars(incvar+6)
   ! ddsen = svars(incvar+7)
   call degradfunc(omegan,domegan,ddomegan,pfn)

   call matmul_uel(nd,dd_copy,ppf,1,1,4)
   pf = ppf(1,1)

   !update crack length
   coordx=0.0_dp
   if (pf.ge.0.95_dp) then
    do i=1,4
       coordx=coordx+nd(1,i)*coords(1,i)
    enddo
    if (coordx.gt.Ac) then
        Ac=coordx
    endif
   endif

   call matmul_uel(nd,ddd,dppf,1,1,4)
   dpf = dppf(1,1)

   call alpha_dd(dalpha,pf,ddalpha,c0)
   call degradfunc(omega,domega,ddomega,pf)

   CALL DF(B,UU,PROPS,Tbar,stiff,PSIP,JELEM,kstep)
   act_psip = psip 
   ! ddse = omega*act_psip - omegan*act_psipn (actual)
   ! ddse = (omega-omegan)*psip

!  Step conditionals 

   if(kstep .eq. 1) then 

   alph=G_MI*PSIP*(1.0_dp-pf)*(1.0_dp-pf)
   ! alph = G_MI*PSIP
   !alph = omega*G_MI*(Hf*SQRT(BETA+1)/G_F + HmI/G_MI + HmII/G_MII)
   !Update fatigue history variable
   if ((alph.ge.alphn).and.(dtime.gt.0.0_dp)) then
    alphB = alphBn+abs(alph-alphn)
    delta_alph = abs(alph-alphn) 
   else
    alphB=alphBn
    delta_alph = 0.0_dp 
   endif      
   
   else

!   The following didn't work
!   alph=G_MI*PSIP*(1.0_dp-pf)*(1.0_dp-pf)
!   ! alphB = alphBn + alph
!   ! alph = G_MI*PSIP
!   if ((alph.ge.alphn).and.(dtime.gt.0.0_dp)) then
!    alphB = alphBn+abs(alph-alphn)
!   else
!    alphB=alphBn
!   endif  

!   The following didn't work
   ! delta_alph = delta_alphn + ddse 

!  Let us try the follows: 
   alphB = alphBn + cycles*delta_alphn !Adding the last delta_alphn
   if(omegan .ne. 0.0_dp) then   
     ratio = (omega/omegan)*(act_psip/act_psipn)  !should be less than or equal to 1.0_dp 
   else 
     ratio = 0.0_dp 
   endif 
   delta_alph = ratio*delta_alphn
   endif

   if (alphB.lt.alphT) then
    Fdeg= 1.0_dp
   else
    Fdeg=(2*alphT/(alphB+alphT))**2
   endif

   if(svars(incvar+1) .gt. psip) then 
      psip = svars(incvar+1)
   else 
      svars(incvar+1) = psip
   endif
   svars(incvar+1) = psip 
   svars(incvar+2) = pf 
   svars(incvar+3) = alph 
   svars(incvar+4) = alphB 
   svars(incvar+5) = delta_alph 
   svars(incvar+6) = act_psip 
   ! svars(incvar+7) = ddse 

   !Storing overall history energy and histories I and II (for visualization)
   UserVar(jelem,1:6,k) = svars(incvar+1:incvar+6)
  
!  Comment one or the other block  
   ph_energy=domega*psip/Fdeg+0.5_dp/lb*dalpha
   d_ph_energy=ddomega*psip/Fdeg+0.5_dp/lb*ddalpha

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call transpose_uel(nd,nd_t,1,4)      
   call transpose_uel(bd,bd_t,2,4)  
   call transpose_uel(b,b_t,3,8)  
!   call matmul_uel(bd_t,bd,bd_t_mul_bd,4,4,2) !NOT USEFUL 
   CALL matmul_uel(BD_T,OMEGA_MAT,BD_T_MUL_OMEGA,4,2,2)
   CALL matmul_uel(BD_T_MUL_OMEGA,BD,BD_T_MUL_BD,4,4,2)  
   call matmul_uel(bd_t_mul_bd,dd,bd_t_mul_bd_mul_dd,4,1,4)
   call matmul_uel(nd_t,nd,nd_t_mul_nd,4,4,1)
   call matmul_uel(b_t,Q_MAT,b_t_mul_d,8,3,3)
   call matmul_uel(b_t_mul_d,b,b_t_mul_d_mul_b,8,8,3)  
                    
   rd1=-w_qp(k)*(ph_energy*nd_t + lb* &   
   bd_t_mul_bd_mul_dd)*det_jacb
   rd=rd+rd1
   
   kd1=w_qp(k)*(d_ph_energy*nd_t_mul_nd+lb* &              
   bd_t_mul_bd)*det_jacb
   kd=kd+kd1

   ku=ku+w_qp(k)* &
   ((omega+0.000001_dp)*b_t_mul_d_mul_b)*det_jacb

!  An Alternative implementation of ru 
   call matmul_uel(b_t_mul_d_mul_b,uu,b_t_mul_d_mul_b_mul_uu,8,1,8)        
   ru = ru - w_qp(k)*((omega+0.000001_dp)*b_t_mul_d_mul_b_mul_uu)*det_jacb

 enddo 

!      call matmul_uel(-ku,uu,ru,8,1,8)
      
 do i=1,4
   do j=1,4
     amatrx(3*i-2,3*j-2)=ku(2*i-1,2*j-1)
     amatrx(3*i-1,3*j-2)=ku(2*i,2*j-1)
     amatrx(3*i-2,3*j-1)=ku(2*i-1,2*j)
     amatrx(3*i-1,3*j-1)=ku(2*i,2*j)
     amatrx(3*i,3*j)=kd(i,j)
   end do
 end do
 do i=1,4
     rhs(3*i-2,1)=ru(2*i-1,1)
     rhs(3*i-1,1)=ru(2*i,1)
     rhs(3*i,1)  =rd(i,1)
 end do

end subroutine element_isotropic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine get_w_qp(w_qp,qpT,qp)
 use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 IMPLICIT NONE 
 !main variables 
 real(dp):: w_qp(4), qpT(2,4), qp(4,2)
 !local variables 
 real(dp):: pos, a_neg 

 w_qp = [1._dp, 1._dp, 1._dp, 1._dp]
 pos = 0.577350269189630_dp 
 a_neg = -0.577350269189630_dp 
 qpT= reshape((/a_neg,a_neg,pos,a_neg, &
     a_neg,pos,pos,pos/), shape(qpT))
 qp = transpose(qpT)
 end subroutine get_w_qp

!==================fem subroutine to provide us the necessary quantities====================== 
 subroutine fem(xi,eta,coords,jelem,k)
 use kvisual 
 use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 implicit none 
 real(dp):: xi,eta,coords(2,4)
 integer:: jelem, k !for integration point not for loops
 !local varibles
 real(dp):: dn_xieta(2,4),jacb(2,2),inv_jacb(2,2) 
 real(dp):: dn_x(4), dn_y(4)
 real(dp):: n(4),det_jacb
 integer:: i
 
 n(1)=1.0_dp/4.0_dp*(1.0_dp-xi)*(1.0_dp-eta)
 n(2)=1.0_dp/4.0_dp*(1.0_dp+xi)*(1.0_dp-eta)
 n(3)=1.0_dp/4.0_dp*(1.0_dp+xi)*(1.0_dp+eta)
 n(4)=1.0_dp/4.0_dp*(1.0_dp-xi)*(1.0_dp+eta)
 ndstor(jelem,k,1:4) = n(1:4) 
 dn_xieta(1,1)=-1.0_dp/4.0_dp*(1.0_dp-eta)
 dn_xieta(1,2)= 1.0_dp/4.0_dp*(1.0_dp-eta)
 dn_xieta(1,3)= 1.0_dp/4.0_dp*(1.0_dp+eta)
 dn_xieta(1,4)=-1.0_dp/4.0_dp*(1.0_dp+eta)
 dn_xieta(2,1)=-1.0_dp/4.0_dp*(1.0_dp-xi)
 dn_xieta(2,2)=-1.0_dp/4.0_dp*(1.0_dp+xi)
 dn_xieta(2,3)= 1.0_dp/4.0_dp*(1.0_dp+xi)
 dn_xieta(2,4)= 1.0_dp/4.0_dp*(1.0_dp-xi)
 jacb=0.0_dp
 do i=1,4
     jacb(1,1)=jacb(1,1)+dn_xieta(1,i)*coords(1,i)
     jacb(1,2)=jacb(1,2)+dn_xieta(1,i)*coords(2,i)
     jacb(2,1)=jacb(2,1)+dn_xieta(2,i)*coords(1,i)
     jacb(2,2)=jacb(2,2)+dn_xieta(2,i)*coords(2,i)
 enddo
 det_jacb=jacb(1,1)*jacb(2,2)-jacb(1,2)*jacb(2,1)
 if(det_jacb .lt. 0.0_dp) then
    write(6,*) 'jacobian is negative of', jelem 
    call xit 
 endif 
 det_jacbstor(jelem,k,1) = det_jacb
 inv_jacb(1,1)=jacb(2,2);  inv_jacb(1,2)=-jacb(1,2)
 inv_jacb(2,1)=-jacb(2,1); inv_jacb(2,2)=jacb(1,1)
 inv_jacb=1.0_dp/det_jacb*inv_jacb
 !     initialize varibles
 ! b=0.0_dp
 ! do i=1,4
 !     nd(1,i)=n(i)
 ! enddo
 do i=1,4
     dn_x(i)=inv_jacb(1,1)*dn_xieta(1,i) &
              +inv_jacb(1,2)*dn_xieta(2,i)
     dn_y(i)=inv_jacb(2,1)*dn_xieta(1,i) &
              +inv_jacb(2,2)*dn_xieta(2,i)
 enddo
 bdstor(jelem,k,1,1:4) = dn_x(1:4)
 bdstor(jelem,k,2,1:4) = dn_y(1:4)

 end subroutine fem 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE CALCULATE_STIFF_COMP(STIFF,PROPS,JELEM)
 use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 IMPLICIT NONE 
 REAL(dp):: STIFF(3,3), PROPS(10)
 REAL(dp):: E_11, E_22, MU_12, NU_12
 INTEGER:: JELEM  
 real(dp) :: nu_21 
 E_11 = PROPS(1) !LONGITUDINAL MODULUS 
 E_22 = PROPS(2) !TRANSVERSE MODULUS
 MU_12 = PROPS(3) !INPLANE SHEAR-MODULUS 
 NU_12 = PROPS(4) !INPLANE'S RATIO
!Calculate nu_21 
nu_21 = NU_12*E_22/E_11
STIFF = 0.0_dp 
stiff(1,1) = E_11/(1.0_dp-nu_12*nu_21) 
stiff(2,2) = E_22/(1.0_dp-nu_12*nu_21)
stiff(1,2) = nu_12*E_22/(1.0_dp-nu_12*nu_21)
stiff(2,1) = nu_12*E_22/(1.0_dp-nu_12*nu_21)
stiff(3,3) = MU_12 
END SUBROUTINE CALCULATE_STIFF_COMP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! subroutine invert_3x3(A, invA)
 ! use iso_fortran_env, only: sp=>real32, dp=>real64
 ! !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 ! implicit none
 ! real(dp), dimension(3,3) :: A, invA
 ! real(dp) :: detA
 ! ! Calculate the determinant of A
 ! detA = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1)+ &
 !        A(1,3)*A(2,1)*A(3,2) - A(1,1)*A(2,3)*A(3,2)- & 
 !        A(1,2)*A(2,1)*A(3,3) - A(1,3)*A(2,2)*A(3,1)
 ! ! Check if the matrix is invertible
 ! if (detA .EQ. 0.0_dp) then
 !   write(6,*) "Matrix is singular, cannot be inverted"
 !   CALL XIT 
 ! endif
 ! ! Calculate the inverse matrix using Cramer's rule
 ! invA(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2)) / detA
 ! invA(1,2) = (A(1,3)*A(3,2) - A(1,2)*A(3,3)) / detA
 ! invA(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2)) / detA
 ! invA(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3)) / detA
 ! invA(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1)) / detA
 ! invA(2,3) = (A(1,3)*A(2,1) - A(1,1)*A(2,3)) / detA
 ! invA(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1)) / detA
 ! invA(3,2) = (A(1,2)*A(3,1) - A(1,1)*A(3,2)) / detA
 ! invA(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1)) / detA
 ! end subroutine invert_3x3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
 SUBROUTINE CALCULATE_TRANS(THETA,JELEM,Tbar,Tbar_trans)
  use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 IMPLICIT NONE 
 REAL(dp):: THETA
 REAL(dp):: C, S
 INTEGER:: JELEM 
 real(dp):: Tbar(3,3), Tbar_trans(3,3)
 real(dp):: a(3,3) 
 C = DCOS(THETA); S = DSIN(THETA)

 !Calculation of T, Tbar and Tbar_trans 
 a = 0.0_dp 
 a(1,1) = C; a(1,2) = S; a(1,3) = 0.0 !First row of rotation matrix 
 a(2,1) = -S; a(2,2) = C; a(2,3) = 0.0 !Second row of rotation matrix 
 a(3,1) = 0.0 ; a(3,2) = 0.0; a(3,3) = 1.0 !Third row of rotation matrix 
 ! T(1,1) = a(1,1)*a(1,1)
 ! T(1,2) = a(1,2)*a(1,2) 
 ! T(1,3) = 2*a(1,1)*a(1,2)
 ! T(2,1) = a(2,1)*a(2,1)
 ! T(2,2) = a(2,2)*a(2,2)
 ! T(2,3) = 2*a(2,1)*a(2,2)
 ! T(3,1) = a(1,1)*a(2,1)
 ! T(3,2) = a(1,2)*a(2,2)
 ! T(3,3) = a(1,1)*a(2,2) + a(1,2)*a(2,1)
 ! CALL invert_3x3(T,T_INV)

 Tbar(1,1) = a(1,1)*a(1,1)
 Tbar(1,2) = a(1,2)*a(1,2) 
 Tbar(1,3) = a(1,1)*a(1,2)
 Tbar(2,1) = a(2,1)*a(2,1)
 Tbar(2,2) = a(2,2)*a(2,2)
 Tbar(2,3) = a(2,1)*a(2,2) 
 Tbar(3,1) = 2*a(1,1)*a(2,1)
 Tbar(3,2) = 2*a(1,2)*a(2,2)
 Tbar(3,3) = a(1,1)*a(2,2) + a(1,2)*a(2,1)  
 call transpose_uel(Tbar,Tbar_trans,3,3)
 END SUBROUTINE CALCULATE_TRANS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE Q_MATRIX(Q_MAT,Tbar_trans,STIFF,Tbar,JELEM) 
 use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 IMPLICIT NONE 
 REAL(dp):: Q_MAT(3,3),Tbar_trans(3,3),STIFF(3,3),Tbar_trans_STIFF(3,3),Tbar(3,3)
 INTEGER:: JELEM 
 CALL matmul_uel(Tbar_trans,STIFF,Tbar_TRANS_STIFF,3,3,3)
 CALL matmul_uel(Tbar_TRANS_STIFF,Tbar,Q_MAT,3,3,3)
 END SUBROUTINE Q_MATRIX

!===============degradation function and its derivations========================
 subroutine degradfunc(omega,domega,ddomega,pf)
 use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 implicit none
 real(dp):: omega, domega, ddomega, pf
 omega = (1.0_dp - pf)*(1.0_dp - pf)
 domega = -2.0_dp*(1.0_dp-pf) 
 ddomega = 2.0_dp 
 end subroutine degradfunc
!=====================r alpha===================================================
 subroutine alpha_dd(dalpha,pf,ddalpha,c0)
 use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 implicit none
 real(dp):: dalpha, pf, ddalpha,c0
 c0 = 2.0_dp 
 dalpha = 2.0_dp*pf 
 ddalpha = 2.0_dp 
 end subroutine alpha_dd  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE DF(B,UU,PROPS,Tbar,stiff,PSIP,JELEM,kstep)
 use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 IMPLICIT NONE 
 REAL(dp):: B(3,8),UU(8,1),PROPS(10),PSIP,Tbar(3,3),stiff(3,3)
 REAL(dp):: EPSI(3,1),EPSI_MAT(3,1),STRESS(3,1),STRESS_MAT(3,1)
 REAL(dp):: PSI_F, PSI_MI, PSI_MII
 REAL(dp):: G_F, G_MI, G_MII, BETA 
 ! real(dp):: Hfn, Hf, HmIn, HmI, HmIIn, HmII
 INTEGER::  JELEM 
 integer:: kstep 
 G_F = PROPS(5) !FIBER FRACTURE TOUGHNESS 
 G_MI = PROPS(6) !MATRIX MODE-I FRACTURE TOUGHNESS 
 G_MII = PROPS(7) !MATRIX MODE-II FRACTURE TOUGHNESS 
 BETA = PROPS(9) !MODEL PENALTY PARAMETER 
 CALL matmul_uel(B,UU,EPSI,3,1,8) !CALCULATE Engineering STRAIN (Global)
 CALL matmul_uel(Tbar, EPSI,EPSI_MAT,3,1,3) !CALCULATE Engineering STRAIN IN MATERIAL CO-ORDINATES
 Call matmul_uel(stiff,EPSI_MAT,stress_mat,3,1,3)
 !Uncomment one of the following for analysis 
 ! Anisotropic decomposition
 ! if(kstep .eq. 1) then  
 ! PSI_F = 0.5*((STRESS_MAT(1,1)+ABS(STRESS_MAT(1,1)))/2)* &
 !              ((EPSI_MAT(1,1)+ABS(EPSI_MAT(1,1)))/2)
 ! PSI_MI = 0.5*((STRESS_MAT(2,1)+ABS(STRESS_MAT(2,1)))/2)* &
 !              ((EPSI_MAT(2,1)+ABS(EPSI_MAT(2,1)))/2)             
 ! PSI_MII = 0.5*STRESS_MAT(3,1)*EPSI_MAT(3,1)
 ! else 
 ! Isotropic decomposition 
 if(kstep .eq. 1) then 
 PSI_F = 0.5*STRESS_MAT(1,1)*EPSI_MAT(1,1)
 PSI_MI = 0.5*STRESS_MAT(2,1)*EPSI_MAT(2,1)          
 PSI_MII = 0.5*STRESS_MAT(3,1)*EPSI_MAT(3,1)
 else 
 PSI_F = 0.5*STRESS_MAT(1,1)*EPSI_MAT(1,1)
 PSI_MI = 0.5*STRESS_MAT(2,1)*EPSI_MAT(2,1)          
 PSI_MII = STRESS_MAT(3,1)*EPSI_MAT(3,1)
 endif 
 ! endif 
 PSIP = PSI_F*SQRT(1._dp+BETA)/G_F + PSI_MI/G_MI + PSI_MII/G_MII
 !PSIP =  PSI_MI/G_MI + PSI_MII/G_MII
 END SUBROUTINE DF

!==================matrix multiplication=======================================
 subroutine matmul_uel(a,b,c,m,n,l)
  use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 implicit none
 integer i, j, k, m, n, l
 real(dp):: a(m,l),b(l,n),c(m,n)
 do i=1,m
   do j=1,n
       c(i,j)=0.0_dp
       do k=1,l
           c(i,j)=c(i,j)+a(i,k)*b(k,j)
       end do
   end do
 end do 
 end subroutine matmul_uel
!==================matrix transpose=======================================
 subroutine transpose_uel(a,a_t,m,n)     
  use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind      
 implicit none
 integer i,j,m,n
 real(dp):: a(m,n),a_t(n,m)
 do i=1,m
     do j=1,n
         a_t(j,i)=a(i,j)
     end do
 end do
 end subroutine transpose_uel
!==================end of the subroutines=======================================

subroutine print_matrix(A, n, m)
use iso_fortran_env, only: sp=>real32, dp=>real64
implicit none
integer, intent(in) :: n, m
real(dp), intent(in), dimension(n, m) :: A
integer :: i, j
do i = 1, n
  write(6, '(25F20.8)') (A(i, j), j = 1, m)
end do
end subroutine print_matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt, &
 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred, &
 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt, &
 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)
 use kvisual
 use iso_fortran_env, only: sp=>real32, dp=>real64
 !use Abaqus_Definitions, only: dp=>abaqus_real_kind
 implicit none 
 character(8) :: cmname
 integer :: ntens, nstatv, nprops, ndi, nshr
 real(dp) :: stress(ntens),statev(nstatv),ddsdde(ntens,ntens), &
 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens), &
 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3), &
 dfgrd0(3,3),dfgrd1(3,3),jstep(4)
 real(dp) :: sse, spd, scd, rpl, drpldt, temp2, dtemp 
 real(dp) :: pnewdt, celent, dtime 
 integer :: noffset, noel, npt, kinc, kspt, layer 
 ddsdde=0.0_dp 
 noffset=noel-nelemp
 statev(1:nstatv-1)=UserVar(noffset,1:nstatv-1,npt)
 statev(nstatv)=Ac
 return
 end subroutine umat 

 