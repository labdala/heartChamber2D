!========================================================================
!
      module mgeometry
!
      implicit none
!
      integer, parameter :: nsd = 2
      integer, parameter :: nf  = 1 ! numero de fontes
!
      integer :: nnodes
      integer :: nelem
      integer :: nedges
      integer :: ndof
!
      integer :: nx
      integer :: ny
!
      integer, parameter :: neng  = 9 ! used for the geometric mapping
      integer, parameter :: nen   = 9 ! number of element nodes (quadrilateral element)
      integer, parameter :: nee   = 4 ! number of element edges (quadrilateral element)
      integer, parameter :: np    = 3 ! 1 para P0 e 3 para P1
!
!     Domain
!
      real(8), dimension(2) :: xint= (/0.d0, 0.04d0/)
      real(8), dimension(2) :: yint= (/0.d0, 0.08d0/)
      real(8) :: h0
!
!     Boundary conditions
!                                       SY   | Pois | correa | scii | zero | degrau
      integer, dimension(6) :: bc= (/   0   ,  0   ,   0    ,   0  ,   0  ,   1   /)
      integer, dimension(4) :: boundc= (/1,1,1,1/) ! 2 p/ tracao e 1 p/ dirichlet
!
!     Pressao prescrita na tracao (caso use boundc 2)? 1-sim, 0-nao
      integer :: ip= 1
!
!     geometria das valvulas (tamanho)
!
      real(8) :: avl,mvl,xmi, xmf, xai, xaf
!
!     types (structures)
!
!     edge
!
      type :: tedge
!
      integer, dimension(4) :: n
      integer, dimension(2) :: e
!
      end type tedge
!
!     quadrilateral
!
      type :: tquad
!
      integer, dimension(nen)   :: n
      integer, dimension(nee)   :: e
!
      real(8), dimension(nee)   :: ngb ! neighbors
!
      integer, dimension(nen)   :: id
!
      end type tquad
!
!     global arrays
!
!     coordinates
!
      real(8), dimension(:,:), allocatable :: coord
!
!     edges
!
      type(tedge), dimension(:), allocatable :: ed
!
!     elements
!
      type(tquad), dimension(:), allocatable :: el
!
      end module mgeometry
!
!----------------------------------------------------------------------
!
      module mgauss
!
      use mgeometry
!
      implicit none
!
      integer, parameter :: nint1dmax=8
      integer, parameter :: nintmax  =125
      real(8), dimension(nintmax,nintmax,nsd) :: xi
      real(8), dimension(nintmax,nintmax)     :: wg
      real(8), dimension(nint1dmax,nint1dmax) :: xig
      real(8), dimension(nint1dmax,nint1dmax) :: wgl
!
      end module
!
!----------------------------------------------------------------------
!
      module mcoeficientes
!
      implicit none
!
      real(8) :: alpha=0.d0
      real(8) :: mu=4.9d-3
      real(8) :: lambda=1.d+10
      real(8) :: rho= 1055.d0
      real(8) :: betapres=1.d-2
      real(8) :: pref=10665.79d+5
      real(8) :: vmin=50e-6
      real(8) :: vmax=120e-6
      real(8) :: tiejec=0.05d0
      real(8) :: tirela=0.3d0
      real(8) :: tifill=0.45d0
      real(8) :: tcicle=1
      real(8) :: t0,t1,t2,t3,t4,t5,t6
      real(8) :: p0=10.d0*133.322365d0
      real(8) :: p1=80.d0*133.322365d0
      real(8) :: p2=120.d0*133.322365d0
      real(8) :: p3=100.d0*133.322365d0
      real(8) :: p4=7.d0*133.322365d0
      real(8) :: p5=5.d0*133.322365d0
      real(8) :: p6=10.d0*133.322365d0
      real(8) :: v0=120e-6
      real(8) :: v1=120e-6
      real(8) :: v2=75e-6
      real(8) :: v3=50e-6
      real(8) :: v4=50e-6
      real(8) :: v5=70e-6
      real(8) :: v6=120e-6
      real(8) :: ava=3.5d-4
      real(8) :: mva=5.d-4
!
      real(8) :: vmed
!
      integer :: ival
!
      end module
!
!----------------------------------------------------------------------
!
      module mtime
!
      implicit none
!
      real(8) :: ti,tf
      integer :: ntot
      real(8) :: dt,tt
!
!     intervalos por trechos
!      
      real(8), dimension(10) :: tic,tfc,dtc
      integer, dimension(10) :: ntc
      integer :: ntrechos
!
      end module
!
!----------------------------------------------------------------------
!
      module mprint
!
      implicit none
!
      real(8) :: tprt_vel,dtprt_vel
      integer :: npvel,nprint_vel
      real(8) :: tprt_pres,dtprt_pres
      integer :: nppres,nprint_pres
      real(8) :: tprt_h,dtprt_h
      integer :: nph,nprint_h
!
      end module
!
!
!----------------------------------------------------------------------
!
      module marrays
!
      implicit none
!
      real(8), dimension(:,:) , allocatable :: p,pa,pk !guarda a solucao das velocidades nos elementos por linhas
      real(8), dimension(:,:) , allocatable :: cpres,cpresa ! guarda os coeficientes (por colunas) que multiplicam as funcoes de base
      real(8), dimension(:,:) , allocatable :: pressao,pressaok,pressaoa !guarda solucao nos nós por linhas
      real(8), dimension(:,:) , allocatable :: thick !guarda valor da espessura nos nós por linhas
      real(8), dimension(:) , allocatable :: thick_aux !guarda valor da espessura nos nós por linhas
!     real(8), dimension(4,8) :: p
!      real(8), dimension(:  ), allocatable :: ph
!
!     lapack: para estrutura de banda
!
      integer :: kl,ku,lda,banda,ld
      real(8), dimension(:,:)  , allocatable :: ab,bb,bk
      real(8), dimension(:,:,:)  , allocatable :: api       ! guarda as matrizes -C1B^T de todos os elementos
      real(8), dimension(:,:,:)  , allocatable :: ag        !guarda C^-1Ge de todos os elementos
!
      real(8) :: pp,ppa,flowr      ! pressao media
!
      end module
!
!----------------------------------------------------------------------
!
      module merro
!
      implicit none
!
      integer :: ind
      integer, parameter :: nexp = 1
      real(8), dimension(nexp,2) :: errol2p, errol2dp, errol2pres
!
      end module
!
!
!
!----------------------------------------------------------------------
!
      module miter
!
      implicit none
!
      integer :: itmax  =100
      real(8) :: epsvel =1.d-2
      real(8) :: epspres=1.d-2
      real(8) :: epsvalvulas=1.d-2
      integer :: indexc
!
      end module
!
!----------------------------------------------------------------------
!
      module maverage
!
      implicit none
!
      real(8) :: vol, presm
      real(8) :: volatot,volmtot
!
      end module
!
!----------------------------------------------------------------------
!
      program navierstokesc
!
!     resolve o problema de navier-stokes compressivel
!
      use mgeometry
      use marrays
      use merro
      use mtime
      use mprint
      use miter
      use mcoeficientes
!
      implicit none

      integer :: n,k, nel, i, ii, jj
      real(8) :: p_out
      real(8) :: calcdt
!
      call setgeometry
!
!     store data for numerical integration
!
      call gausstable
!
      nx= 20
      ny= 20
!
      p_out = 0.d0!pref/2.d0
      ival= 0
!
      call settime
!
!     inicializa os contadores da impressao
!
      call initprt(ti,tf)
!
!     read mesh information
!
      call genmesh
!
!     define ndof:
!
      ndof = 2*nnodes !2*nelem + 1
!
!     para lapack
!
      call calcband
!
      ku = banda
      kl = ku
!
      lda = ku+1+kl+kl
      ld  = kl+ku+1  ! posicao da diagonal
!
      allocate(ab(lda,ndof))
      allocate(bb(ndof,1))
      allocate(bk(ndof,1))
      allocate(api(nelem,np,2*nen))
      allocate(ag(nelem,np,1))
      ab  = 0.d0
      bb  = 0.d0
      api = 0.d0
      ag  = 0.d0
!
!     aloca os arrays globais
!
      allocate(p (nelem,2*nen))
      allocate(pa(nelem,2*nen))
      allocate(pk(nelem,2*nen))
      allocate(cpres (np,nelem))
      allocate(cpresa(np,nelem))
      allocate(pressao (nelem,nen))
      allocate(pressaok(nelem,nen))
      allocate(pressaoa(nelem,nen))
      allocate(thick(nelem,nen))
      allocate(thick_aux(4*nelem))
!
!     gera o .case para imprimir no paraview
!
      call geracase_vel (npvel )
      call geracase_pres(nppres)
      call geracase_h   (nph  )

      call paraview_geom
!
!     inicializa o tempo
!
      tt=ti
!
!     inicializa a solucao
!
      p      =0.d0
      pressao=1333.22d0
      cpres  =0.d0
      cpresa =0.d0
      thick  =0.d0 
      thick_aux=0.d0
!
!     imprime as condicoes iniciais
!
      call paraview_vel(nprint_vel)
      tprt_vel=tprt_vel+dtprt_vel
!
      call paraview_pres(nprint_pres)
      tprt_pres=tprt_pres+dtprt_pres
!
!      call paraview_h(nprint_h)
!      tprt_h=tprt_h+dtprt_h
!
!      bb=0.d0
!
!      do n=1,ntot
      n=1
      do while(tt.le.tf)
!
!     atualiza o tempo
!
      dt=calcdt(tt)
!
      tt=tt+dt
!
!     armazena a solucao anterior
!
      pa=p
      bb=0.d0
      cpresa=cpres
      pressaoa=pressao
!
      do jj=ny,ny
      do ii=1 ,nx
      nel=ii+nx*(jj-1)
!
!     Aortica: Compressao isovol + Ejecao
!
      if(tt.le.t1)then
!
      if(coord(el(nel)%n(7),1).ge.xai.and.coord(el(nel)%n(7),1).le.xaf) then
      el(nel)%id(3) =6
      el(nel)%id(4) =6
      el(nel)%id(7) =6
      endif
!
      elseif(tt.gt.t1.and.tt.le.t3)then
!
      if(coord(el(nel)%n(7),1).ge.xai.and.coord(el(nel)%n(7),1).le.xaf) then
      el(nel)%id(3) =6
      el(nel)%id(4) =6
      el(nel)%id(7) =6
      endif
!
!     Mitral: Relaxamento isovol
!
      elseif(tt.gt.t3.and.tt.le.t4) then
!
      if(coord(el(nel)%n(7),1).ge.xai.and.coord(el(nel)%n(7),1).le.xaf) then
      el(nel)%id(3) =1
      el(nel)%id(4) =1
      el(nel)%id(7) =1
!
      elseif(coord(el(nel)%n(7),1).ge.xmi.and.coord(el(nel)%n(7),1).le.xmf) then
      el(nel)%id(3) =6
      el(nel)%id(4) =6
      el(nel)%id(7) =6
!
      end if
!
!     Mitral: enchimento
!
      elseif(tt.gt.t4) then
!
      if(coord(el(nel)%n(7),1).ge.xai.and.coord(el(nel)%n(7),1).le.xaf) then
      el(nel)%id(3) =1
      el(nel)%id(4) =1
      el(nel)%id(7) =1
!
      elseif(coord(el(nel)%n(7),1).ge.xmi.and.coord(el(nel)%n(7),1).le.xmf) then
      el(nel)%id(3) =6
      el(nel)%id(4) =6
      el(nel)%id(7) =6
!
      endif
!
      else
      write(*,*) "tempo maior que um ciclo 6"
      stop
!
      endif
!
      end do ! ii
      end do ! jj
!
      if(tt.gt.t1.and.tt.le.t3.or.tt.gt.t4.and.tt.le.t6) then
!
!     inicializa o processo iterativo sequencial
!
      do k=1,itmax
!      write(*,*) k
!
!      if(n.eq.1) cpresa=cpres
!
!     inicializa os arrays
!
      bk=bb
      pk=p
      pressaok=pressao
!
      ab=0.d0
      bb=0.d0
!
      call assemble
!
!     solving the global system
!
      call solve
!
!     postprocessing the solution
!
      call postprocessing
!
!     checa a convergencia
!
      call convcheck
!
      if(indexc.eq.1) then
      write(*,"('n=',i5,' tn=',f12.5,' convergiu com',i5,' iteracoes')") n,tt,k
      exit
      end if
!
      end do ! k
!
      if(k.ge.itmax) then
      write(*,*) " ATENCAO: numero maximo de iteracoes atingido: ",k
   !   stop
      end if
!
!     caso contrario, esta no trecho isovolumetrico
!
      else
!
      call solprescrita
      write(*,"('n=',i5,' tn=',f12.5,' trecho isovolumetrico')") n,tt
!
      end if
!
      open(unit=63,file="volpres.dat",status="unknown")
      open(unit=64,file="vvalvulas.dat",status="unknown")
!
!     calcula e impime as variaveis medias
!
!      write(*,*) tt
!
!     imprime a pressao media e o volume da camara
!
      call volpres
!
!     calcula os fluxos nas valvulas
!
      call vvalvulas
!
!.... imprime a velocidade
!
      if(tt.ge.tprt_vel) then
!
!     calculando a velocidade no centro
!
!      call vel_c(nx,ny,nz,numel,nsd,nf,v,vc)
!
      call paraview_vel(nprint_vel)
      tprt_vel=tprt_vel+dtprt_vel
      end if
!      call prtgnuplot
!
!     imprime a espessura
!
      if(tt.ge.tprt_h) then
      call paraview_h(nprint_h)
      tprt_h=tprt_h+dtprt_h
      end if
!
!     imprime a pressao
!
      if(tt.ge.tprt_pres) then
      call paraview_pres(nprint_pres)
      tprt_pres=tprt_pres+dtprt_pres
      end if
!
!
!     calcula a norma l2
!
!      if(i.eq.1) then
!      write(*,*) 'norma l2'
!      end if
!      call norms
!
!      end do ! n
      n=n+1
      end do ! while
!
!     postprocessing the solution
!
      call postprocessing
!
!      call volpres
      call vvalvulas
!
!     checa o balanco
!
      call vcheck
!
      close(unit=63)
      close(unit=64)
!
!     imprime a solucao no instante final
!
      call prtgnuplot
!
      deallocate(coord,el,ed,ab,bb,bk,p,pa,pk,api, ag, cpres,cpresa, pressao,pressaok,thick, thick_aux)
!
      end program
!
!----------------------------------------------------------------------
!
      subroutine solprescrita
!
      use mgeometry
      use marrays
!
      implicit none
!
      real(8), dimension(nsd)     :: xx
      real(8) :: pr
!
      xx=0.d0
!
!     velocidade
!
      p = 0.d0
!
!     pressao
!
      pressao=pr(xx)
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine settime
!
      use mcoeficientes
      use mtime

      implicit none
!
      real(8) :: eps
      integer :: i
!
      ntrechos=6
!
!
      eps=(tf-ti)/1.d4
!
      t0=0.d0
      t1=0.05d0+eps
      t2=0.2d0+eps
      t3=0.3d0+eps
      t4=0.45d0+eps
      t5=0.6d0+eps
      t6=1.d0
!
      ti=0.00d-0
      tf=0.99d0!t6!0.5d0!t2!0.2d0!t6!0.0d-0
!
!      ntot=1000!210
!      dt=(tf-ti)/ntot
!
!      por trechos
!      
      tic(1)=ti
      tfc(1)=t1
      ntc(1)=10
!      
      tic(2)=tfc(1)
      tfc(2)=t2
      ntc(2)=20
!      
      tic(3)=tfc(2)
      tfc(3)=t3
      ntc(3)=20
!      
      tic(4)=tfc(3)
      tfc(4)=t4
      ntc(4)=20
!      
      tic(5)=tfc(4)
      tfc(5)=t5
      ntc(5)=25
!      
      tic(6)=tfc(5)
      tfc(6)=tf
      ntc(6)=25
!      
      do i=1,ntrechos
      dtc(i)=(tfc(i)-tic(i))/ntc(i)
      end do ! i
!
      end subroutine
!
!----------------------------------------------------------------------
!
      function calcdt(t)
!
      use mtime
      implicit none
!
      real(8) :: t,calcdt
      integer :: i
!
      do i=1,ntrechos
      if(t.ge.tic(i).and.t.lt.tfc(i)) then
      calcdt=dtc(i)
      exit
      end if
      end do ! i      
!      
      if(i.gt.ntrechos) then
      write(*,*) "Tempo fora dos intervalos em calcdt"
      stop
      end if
!      
      end function
!
!----------------------------------------------------------------------
!
      subroutine setgeometry
!
      use mgeometry
      use mcoeficientes

      implicit none
!
      real(8) :: xi, xl, yi, yl,lx,ly
      real(8) :: eps=1.d-8
!
      xi = xint(1)
      yi = yint(1)
      xl = xint(2)
      yl = yint(2)
!
      ly   = yint(2)-yint(1)
      lx   = xint(2)-xint(1)
      vmed = (vmin+vmax)/2.d0
      h0   = vmed/(lx*ly)
!
      avl=ava/h0
      mvl=mva/h0
!
!     limites da aortica
!
      xai=xint(1)
      xaf=xint(1)+avl
!
!     limites da mitral
!
      xmi=xint(2)-mvl
      xmf=xint(2)
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine convcheck
!
!     checa a convergencia do sistema da hidrodinamica usando o multiplicador
!
      use miter
      use mgeometry
!
      use mgauss
      use marrays
      use mcoeficientes
!

      implicit none
!
      real(8) :: enorm,bnorm,aux
      integer :: i,k
      real(8) :: eflux,flux,fluxk,efluxe,ff,ffe
      integer :: nint=64,nint1d=8
      real(8), dimension(nen,nsd) :: xe
      real(8), dimension(4)       :: ee
      integer :: ii,jj,j,l,n1,n2,nel
      real(8) :: xil1d,phi2d
      real(8) :: f4,f3,f7
      real(8), dimension(nsd)     :: xil,xxl
!
!
!     indice de convergencia
!
!     indexc=0 = nao convergiu
!     indexc=1 = convergiu
!
      indexc=0
!
      enorm=0.d0
      bnorm=0.d0
      do i=1,ndof
      enorm=enorm+(bb(i,1)-bk(i,1))**2.d0
      bnorm=bnorm+(bb(i,1))**2.d0
      end do  ! i
      enorm=dsqrt(enorm)
      bnorm=dsqrt(bnorm)
!      aux=min(bnorm,1.d0)
      aux=bnorm ! (relativo)
!      aux =1d0 ! (absoluto)
!
!     testa o fluxo na aortica e na mitral
!
      eflux=0.d0
      ff   =0.d0
!      
      do jj=ny,ny
      do ii=1 ,nx
      nel=ii+nx*(jj-1)
!
      do j=1,nen
      do k=1,nsd
      xe(j,k) = coord(el(nel)%n(j),k)
      end do ! k
      end do ! j

      ee(1)=dsqrt((xe(1,1)-xe(4,1))**2.d0+(xe(1,2)-xe(4,2))**2.d0)
      ee(2)=dsqrt((xe(2,1)-xe(1,1))**2.d0+(xe(2,2)-xe(1,2))**2.d0)
      ee(3)=dsqrt((xe(3,1)-xe(2,1))**2.d0+(xe(3,2)-xe(2,2))**2.d0)
      ee(4)=dsqrt((xe(4,1)-xe(3,1))**2.d0+(xe(4,2)-xe(3,2))**2.d0)
!
!     erro no fluxo da aortica e mitral, entre iteracoes
!
      if(coord(el(nel)%n(7),1).ge.xai.and.coord(el(nel)%n(7),1).le.xaf) then
!
      efluxe=0.d0
      ffe   =0.d0
!
      do l=1,nint1d
!
      xil1d = xig(l,nint1d)
      xil(1) = xil1d
      xil(2) = 1.d0
!
!     Real coordinates of the integration point
!
      xxl = 0.d0
!
      do n1=1,neng
      do n2=1,nsd
      xxl(n2) = xxl(n2) + phi2d(n1,xil,neng)*xe(n1,n2)
      end do ! n1
      end do ! n2
!
!     fluxo normal no ponto de integracao
!
      f4=pk(nel, 8)
      f7=pk(nel,14)
      f3=pk(nel, 6)
      flux =f4*phi2d(4,xil,nen)+f7*phi2d(7,xil,nen)+f3*phi2d(3,xil,nen)
!
      f4=p(nel, 8)
      f7=p(nel,14)
      f3=p(nel, 6)
      fluxk=f4*phi2d(4,xil,nen)+f7*phi2d(7,xil,nen)+f3*phi2d(3,xil,nen)
!
      efluxe=efluxe+((flux-fluxk)**2.d0)*ee(4)/2.d0*wgl(l,nint1d)
      ffe   =ffe+flux**2.d0*ee(4)/2.d0*wgl(l,nint1d)
!
      end do ! l
!
      eflux=eflux+efluxe
      ff   =ff   +ffe
!
      end if
!
!      mitral
!
      if(coord(el(nel)%n(7),1).ge.xmi.and.coord(el(nel)%n(7),1).le.xmf) then
!
      efluxe=0.d0
      ffe   =0.d0
!
      do l=1,nint1d
!
      xil1d = xig(l,nint1d)
      xil(1) = xil1d
      xil(2) = 1.d0
!
!     Real coordinates of the integration point
!
      xxl = 0.d0
!
      do n1=1,neng
      do n2=1,nsd
      xxl(n2) = xxl(n2) + phi2d(n1,xil,neng)*xe(n1,n2)
      end do ! n1
      end do ! n2
!
!     fluxo normal no ponto de integracao
!
      f4=pk(nel, 8)
      f7=pk(nel,14)
      f3=pk(nel, 6)
      flux =f4*phi2d(4,xil,nen)+f7*phi2d(7,xil,nen)+f3*phi2d(3,xil,nen)
!
      f4=p(nel, 8)
      f7=p(nel,14)
      f3=p(nel, 6)
      fluxk=f4*phi2d(4,xil,nen)+f7*phi2d(7,xil,nen)+f3*phi2d(3,xil,nen)
!
      efluxe=efluxe+((flux-fluxk)**2.d0)*ee(4)/2.d0*wgl(l,nint1d)
      ffe   =ffe+flux**2.d0*ee(4)/2.d0*wgl(l,nint1d)
!
      end do ! l
!
      eflux=eflux+efluxe
      ff   =ff   +ffe
!
      end if
      end do ! ii
      end do ! jj
!
!     erro do fluxo na aortica e na mitral
!
      eflux=dsqrt(eflux)
      ff   =dsqrt(ff)
!
      write(*,*) 'aortica e mitral',eflux,epsvalvulas*ff
!
!     testa 
!
      if(enorm.le.epsvel*aux.and.eflux.le.epsvalvulas*ff) then
      indexc=1
      else
      indexc=0
      end if 
!      
      end
!
!
!----------------------------------------------------------------------
!
      subroutine geracase_pres(ns)
      implicit none
!
      integer :: ns
      integer :: i,nz=0,nu=1
!
      character(len=5) :: name1
      character(len=128) :: name2,sol
      character(len=8)   :: c
      real :: x=1
!
      open(unit=35,file="./paraview/pres/pres.case",status="unknown")
!
      write(35,'(a)') 'FORMAT'
      write(35,'(a)') 'type: ensight'
      write(35,'(a)')
      write(35,'(a)') 'GEOMETRY'
      write(35,'(a)') 'model: ../geo/solucao.geo'
      write(35,'(a)')
      write(35,'(a)') 'VARIABLE'
!
      name1="P"
!
      name2="./pres.***"
      write(35,'(1000a)') 'scalar per element: ',name1,name2
!
      write(35,'(a)')
      write(35,'(a)') 'TIME'
      write(35,'(a,i1)') 'time set: ',nu
      write(35,'(a,i3)') 'number of steps: ',ns
      write(35,'(a,i1)') 'filename start number: ',nz
      write(35,'(a,i1)') 'filename increment: ',nu
      write(35,'(a,i5)') 'time values: '
      do i=1,ns
      write(35,'(f10.5)') 1.0*(i-1)
      end do
!
      close(35)
!
      end subroutine

!----------------------------------------------------------------------
!
      subroutine geracase_h(ns)
      implicit none
!
      integer :: ns
      integer :: i,nz=0,nu=1
!
      character(len=5) :: name1
      character(len=128) :: name2,sol
      character(len=8)   :: c
      real :: x=1
!
      open(unit=35,file="./paraview/h/h.case",status="unknown")
!
      write(35,'(a)') 'FORMAT'
      write(35,'(a)') 'type: ensight'
      write(35,'(a)')
      write(35,'(a)') 'GEOMETRY'
      write(35,'(a)') 'model: ../geo/solucao.geo'
      write(35,'(a)')
      write(35,'(a)') 'VARIABLE'
!
      name1="H"
!
      name2="./h.***"
      write(35,'(1000a)') 'scalar per element: ',name1,name2
!
      write(35,'(a)')
      write(35,'(a)') 'TIME'
      write(35,'(a,i1)') 'time set: ',nu
      write(35,'(a,i3)') 'number of steps: ',ns
      write(35,'(a,i1)') 'filename start number: ',nz
      write(35,'(a,i1)') 'filename increment: ',nu
      write(35,'(a,i5)') 'time values: '
      do i=1,ns
      write(35,'(f10.5)') 1.0*(i-1)
      end do
!
      close(35)
!
      end subroutine

!----------------------------------------------------------------------
!
!
      subroutine solve
!
!     solve the global system AX=B
!
      use mgeometry
      use marrays
!
      implicit none
!
      integer :: i,info,ipiv(ndof)
!
!     subrotinas lapack
!
!     xyyzzz
!
!     x   = d   (double)
!     yy  = gb  (general band)
!     zzz = trf (factorize)
!     zzz = trs (substitution)
!
      call dgbsv(ndof,kl,ku,1,ab,lda,ipiv,bb,ndof,info)
!
!     APAGUEI VARIOS COMENTARIOS DO ORIGINAL, pois nao devia!
!
      end subroutine solve
!
!      
!----------------------------------------------------------------------
!
      subroutine paraview_pres(n)
      use marrays
      use mgeometry

      implicit none
!
!     imprime solucoes escalares
!
      character(len=128) :: name,sol
      character(len=8)   :: c
      integer :: i,n,ii
      real(4) ::x
      real(4), dimension(4*nelem) :: paux
!
      x=0.001
!
      sol="./paraview/pres/pres"
!
      write(c,"(f4.3)") x*n
      c=adjustl(c)
      name=trim(sol)//c !trim(c)
!    name=adjustl(trim(name))
!
      open(unit=45,file=name,status="unknown")
!
      write(45,"('Ensight Scalar passo ',i5)") n
      write(45,"('part 1')")
      write(45,"('hexa8')")
!
!   imprime as coordenadas
!
      paux=0.0
      do i=1,nelem
      ii=(i-1)*4
!
!     primeiro subelemento
!
      paux(ii+1)=real(0.25d0*((pressao(i, 1)+pressao(i, 5)+pressao(i, 9)+pressao(i, 8))))
!
!     segundo subelemento
!
      paux(ii+2)=real(0.25d0*((pressao(i, 5)+pressao(i, 2)+pressao(i, 6)+pressao(i, 9))))
!
!     terceiro subelemento
!
      paux(ii+3)=real(0.25d0*((pressao(i, 8)+pressao(i, 9)+pressao(i, 7)+pressao(i, 4))))
!
!     quarto subelemento
!
      paux(ii+4)=real(0.25d0*((pressao(i, 9)+pressao(i, 6)+pressao(i, 3)+pressao(i, 7))))
!
      end do ! i

      write(45,"(6(e12.5))") (real(paux(i)),i=1,4*nelem) !*thick_aux(i)/0.036d0
!
      close(45)
!
      n=n+1
!
      end subroutine
!----------------------------------------------------------------------
!
      subroutine paraview_h(n)
      use marrays
      use mgeometry

      implicit none
!
!     imprime solucoes escalares
!
      character(len=128) :: name,sol
      character(len=8)   :: c
      integer :: i,n,ii
      real(4) ::x
!      real(4), dimension(4*nelem) :: thick_aux !cada elemento da malha eh dividido em 4
!
      x=0.001
!
      sol="./paraview/h/h"
!
      write(c,"(f4.3)") x*n
      c=adjustl(c)
      name=trim(sol)//c !trim(c)
!    name=adjustl(trim(name))
!
      open(unit=45,file=name,status="unknown")
!
      write(45,"('Ensight Scalar passo ',i5)") n
      write(45,"('part 1')")
      write(45,"('hexa8')")
!
!   imprime as coordenadas
!
      thick_aux=0.0
      do i=1,nelem
      ii=(i-1)*4
!
!     primeiro subelemento
!
      thick_aux(ii+1)=real(0.25d0*((thick(i, 1)+thick(i, 5)+thick(i, 9)+thick(i, 8))))
!
!     segundo subelemento
!
      thick_aux(ii+2)=real(0.25d0*((thick(i, 5)+thick(i, 2)+thick(i, 6)+thick(i, 9))))
!
!     terceiro subelemento
!
      thick_aux(ii+3)=real(0.25d0*((thick(i, 8)+thick(i, 9)+thick(i, 7)+thick(i, 4))))
!
!     quarto subelemento
!
      thick_aux(ii+4)=real(0.25d0*((thick(i, 9)+thick(i, 6)+thick(i, 3)+thick(i, 7))))
!
      end do ! i

      write(45,"(6(e12.5))") (real(thick_aux(i)),i=1,4*nelem)
!
      close(45)
!
      n=n+1
!
      end subroutine
!----------------------------------------------------------------------
!
      subroutine paraview_vel(n)
      use marrays
      use mgeometry
      implicit none
!
      integer :: i,n
      character(len=128) :: name,sol
      character(len=8)   :: c
      real(4) ::x
      real(8) :: vc1,vc2,vc3
      real(8) :: vc1a,vc2a,vc3a
!
      x=0.001
      vc3 =0.d0
      vc3a=0.d0
!
      sol="./paraview/vel/vel"
!
      write(c,"(f4.3)") x*n
      c=adjustl(c)
      name=trim(sol)//c !trim(c)
!    name=adjustl(trim(name))
!
      open(unit=45,file=name,status="unknown")
!
      write(45,"('Ensight Vector passo ',i5)") n
      write(45,"('part 1')")
      write(45,"('hexa8')")
!
!     imprime
!
!      write(45,"(6(e12.5))") (real(p(i,17)),real(p(i,18)),x,i=1,nelem)

      do i=1,nelem
!
!     primeiro subelemento
!
      vc1=real(0.25d0*((p(i, 1)+p(i, 9)+p(i,17)+p(i,15))))
      vc2=real(0.25d0*((p(i, 2)+p(i,10)+p(i,18)+p(i,16))))
!
!     segundo subelemento
!
      vc1a=real(0.25d0*((p(i, 9)+p(i, 3)+p(i,11)+p(i,17))))
      vc2a=real(0.25d0*((p(i,10)+p(i, 4)+p(i,12)+p(i,18))))
!
      write(45,"(6(e12.5))") vc1,vc2,vc3,vc1a,vc2a,vc3a
!
!     terceiro subelemento
!
      vc1=real(0.25d0*((p(i,15)+p(i,17)+p(i,13)+p(i, 7))))
      vc2=real(0.25d0*((p(i,16)+p(i,18)+p(i,14)+p(i, 8))))
!
!     quarto subelemento
!
      vc1a=real(0.25d0*((p(i,17)+p(i,11)+p(i, 5)+p(i,13))))
      vc2a=real(0.25d0*((p(i,18)+p(i,12)+p(i, 6)+p(i,14))))
!
      write(45,"(6(e12.5))") vc1,vc2,vc3,vc1a,vc2a,vc3a

      end do !i
!
      close(45)
!
      n=n+1
!
      end subroutine
!

!----------------------------------------------------------------------
!
      subroutine geracase_vel(ns)
      implicit none
!
      integer :: ns
      integer :: i,nz=0,nu=1
!
      character(len=5) :: name1
      character(len=128) :: name2,sol
      character(len=8)   :: c
      real :: x=1
!
      open(unit=35,file="./paraview/vel/vel.case",status="unknown")
!
      write(35,'(a)') 'FORMAT'
      write(35,'(a)') 'type: ensight'
      write(35,'(a)')
      write(35,'(a)') 'GEOMETRY'
      write(35,'(a)') 'model: ../geo/solucao.geo'
      write(35,'(a)')
      write(35,'(a)') 'VARIABLE'
!
      name1="V"
!
      name2="./vel.***"
      write(35,'(1000a)') 'vector per element: ',name1,name2
!
      write(35,'(a)')
      write(35,'(a)') 'TIME'
      write(35,'(a,i1)') 'time set: ',nu
      write(35,'(a,i3)') 'number of steps: ',ns
      write(35,'(a,i1)') 'filename start number: ',nz
      write(35,'(a,i1)') 'filename increment: ',nu
      write(35,'(a,i5)') 'time values: '
      do i=1,ns
      write(35,'(f10.5)') 1.0*(i-1)
      end do
!
      close(35)
!
      end subroutine

!
!----------------------------------------------------------------------
!
      subroutine paraview_geom

      use marrays
      use mgeometry
      implicit none
!
      real(8) :: zz
	  integer :: i,j
!
      open(unit=35,file="./paraview/geo/solucao.geo",status="unknown")

	  write(35,'(a)')'Title1'
      write(35,'(a)')'Title2'
      write(35,'(a)')'node id given'
      write(35,'(a)')'element id given'
      write(35,'(a)')'coordinates'
      write(35,'(i8)') 2*nnodes

!
!     imprime as coordenadas
!
      zz=0.d0
      do i=1,nnodes
      write(35,"(i8,3(e12.5))") i,coord(i,1),coord(i,2),zz
      end do ! i
!
!     imprime as coordenadas
!
      zz=0.001d0
      do i=1,nnodes
      write(35,"(i8,3(e12.5))") i+nnodes,coord(i,1),coord(i,2),zz
      end do ! i
!

	  write(35,'(A,/,A,/,A,/,I8)') 'part 1'           , &
                                              'malha'            , &
                                              'hexa8'            ,4*nelem
!
!     conectividades
!
!      do i=1,nelem
!      write(35,"(9i8)") i,el(i)%n(1),el(i)%n(2),el(i)%n(3),el(i)%n(4)
!      end do ! i
!
      j=1
      do i=1,nelem
!
!     primeiro subelemento
!
      write(35,"(9i8)") j,el(i)%n(1),el(i)%n(5),el(i)%n(9),el(i)%n(8), &
                          el(i)%n(1)+nnodes,el(i)%n(5)+nnodes,el(i)%n(9)+nnodes,el(i)%n(8)+nnodes
      j=j+1
!
!     segundo subelemento
!
      write(35,"(9i8)") j,el(i)%n(5),el(i)%n(2),el(i)%n(6),el(i)%n(9), &
                          el(i)%n(5)+nnodes,el(i)%n(2)+nnodes,el(i)%n(6)+nnodes,el(i)%n(9)+nnodes
      j=j+1
!
!     terceiro subelemento
!
      write(35,"(9i8)") j,el(i)%n(8),el(i)%n(9),el(i)%n(7),el(i)%n(4), &
                          el(i)%n(8)+nnodes,el(i)%n(9)+nnodes,el(i)%n(7)+nnodes,el(i)%n(4)+nnodes
      j=j+1
!
!     quarto subelemento
!
      write(35,"(9i8)") j,el(i)%n(9),el(i)%n(6),el(i)%n(3),el(i)%n(7), &
                          el(i)%n(9)+nnodes,el(i)%n(6)+nnodes,el(i)%n(3)+nnodes,el(i)%n(7)+nnodes
      j=j+1
!
      end do ! i
      close(35)
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine initprt(ti,tf)
!
      use mprint
!
      implicit none
      real(8) :: ti,tf,dt
!
      dt=tf-ti
!
!     numero de impressoes
!
      npvel   = 50
      nppres  = 50
      nph    = 50
!
!     contadores
!
      nprint_vel   =0
      nprint_pres  =0
      nprint_h     =0
!
!.... tamanhos dos intervalos de impressoes
!
      dtprt_vel   = dt/npvel
      dtprt_pres  = dt/nppres
      dtprt_h     = dt/nph
!
!.... inicializa os contadores de impressao
!
      tprt_vel   = 0.d0 ! dtprt_vel
      tprt_pres  = 0.d0
      tprt_h     = 0.d0
!
      end
!
!----------------------------------------------------------------------
!
      subroutine calcband
!
      use mgeometry
      use marrays
      implicit none
!
      select case(nen)
!
!      case(1)
!
!      banda=1
!
      case(4)
!
      banda=2*nx+5
!
      case(9)
!
      banda=8*nx+9
!
!      case(8)
!
!      banda=2*(2*nx+1)+2
!
!      case(16)
!
!      banda=3*(3*nx+1)+3
!
      case default
!
      write(*,*) 'element degree not defined'
      stop
!
      end select
!
      end subroutine
!----------------------------------------------------------------------
!
      function ibanda(i,j)
!
      use marrays
!
      implicit none
!
      integer             :: ibanda
      integer, intent(in) :: i,j
!
      ibanda = kl+ku+1+i-j
!
      end function
!
!----------------------------------------------------------------------
!
      function jbanda(i,j)
!
      use marrays
!
      implicit none
!
      integer             :: jbanda
      integer, intent(in) :: i,j
!
      jbanda = j
!
      end function
!
!----------------------------------------------------------------------
!
      subroutine genmesh
!
      use mgeometry
      use mcoeficientes
      use mtime
      implicit none
!
      integer :: i,j,l,k,l1,l2,l3,l4
      real(8) :: ly, lx
      real(8) :: hx,hy,xl,yl,xi,yi, um4, um2
!
      xi = xint(1)
      yi = yint(1)
      xl = xint(2)
      yl = yint(2)
!
      select case(nen)
!
      case(4)
!
!     bilinear
!
      hx = (xl-xi)/nx
      hy = (yl-yi)/ny
!
      nnodes = (nx+1)*(ny+1)
      nelem  = nx*ny
      nedges = 2*nx*ny+nx+ny
!
!     allocating arrays
!
      allocate(coord(nnodes,nsd))
      allocate(el(nelem))
      allocate(ed(nedges))
!
!     coordinates
!
      do j=1,ny+1
!
      do i=1,nx+1
!
      l=i+(nx+1)*(j-1)
      coord( l,1) = xi+(i-1)*hx
      coord( l,2) = yi+(j-1)*hy
!
      end do ! i
!
      end do ! j
!
!     conectivities: nodes
!
      do j=1,ny
!
      do i=1,nx
!
      l=i+(nx)*(j-1)
      el(l)%n(1)=l+j-1
      el(l)%n(2)=l+j
      el(l)%n(3)=l+j+nx+1
      el(l)%n(4)=l+j+nx
!
      do k = 1, nen
!
      if (coord(el(l)%n(k),1)==xi) then
      el(l)%id(k) = boundc(1)
      else if (coord(el(l)%n(k),1)==xl) then
      el(l)%id(k) = boundc(3)
      else if (coord(el(l)%n(k),2)==yi) then
      el(l)%id(k) = boundc(2)
      else if (coord(el(l)%n(k),2)==yl) then
      el(l)%id(k) = boundc(4)
      else
      el(l)%id(k) = 0
      end if
!

!     Impor que nos do canto da malha sejam sempre de Dirichlet nos encontros 2 e 1
!
      if (coord(el(l)%n(k),1)==xi .and. coord(el(l)%n(k),2)==yi) then
      if (boundc(1)==1 .or. boundc(2)==1) el(l)%id(k) = 1
!
      else if(coord(el(l)%n(k),1)==xl .and. coord(el(l)%n(k),2)==yi) then
      if (boundc(2)==1 .or. boundc(3)==1) el(l)%id(k) = 1
!
      else if(coord(el(l)%n(k),1)==xl .and. coord(el(l)%n(k),2)==yl) then
      if (boundc(3)==1 .or. boundc(4)==1) el(l)%id(k) = 1
!
      else if(coord(el(l)%n(k),1)==xi .and. coord(el(l)%n(k),2)==yl) then
      if (boundc(4)==1 .or. boundc(1)==1) el(l)%id(k) = 1
!
      end if
!
!      print *,  el(l)%id(k)
      end do ! k
!
      end do ! i
!
      end do ! j
!
!      do i=1,nelem
!      write(*,*) (el(i)%n(j),j=1,4)
!      end do
!      stop
!
!........................................
!
      case(9)
!
!     biquadratic

      hx=(xl-xi)/nx/2.d0
      hy=(yl-yi)/ny/2.d0
!
      nnodes= (2*nx+1)*(2*ny+1)
      nelem = nx*ny
!      nedges= 2*nx*ny+nx+ny
!
!     allocating arrays
!
      allocate(coord(nnodes,nsd))
      allocate(el(nelem))
      allocate(ed(nedges))
!
!     coordinates
!
      do j=1,(2*ny+1)
!
      do i=1,(2*nx+1)
!
      l=i+(2*nx+1)*(j-1)
      coord( l,1) = xi+(i-1)*hx
      coord( l,2) = yi+(j-1)*hy
!
      end do ! i
!
      end do ! j
!
!      do i=1,nnodes
!      write(*,*) (coord(i,j),j=1,2)
!      end do
!      stop
!
!     conectivities: nodes
!
      do j=1,ny
!
      do i=1,nx
!
      l=i+(nx)*(j-1)
      el(l)%n(1)=2*(j-1)*(2*nx+1)+2*i-1
      el(l)%n(5)=el(l)%n(1)+1
      el(l)%n(2)=el(l)%n(1)+2
!
      el(l)%n(8)=el(l)%n(1)+(2*nx+1)
      el(l)%n(9)=el(l)%n(8)+1
      el(l)%n(6)=el(l)%n(8)+2
!
      el(l)%n(4)=el(l)%n(1)+2*(2*nx+1)
      el(l)%n(7)=el(l)%n(4)+1
      el(l)%n(3)=el(l)%n(4)+2
!
      ! contorno
!
     do k = 1, nen
!
      if (coord(el(l)%n(k),1)==xi) then
      el(l)%id(k) = boundc(1)
      else if (coord(el(l)%n(k),1)==xl) then
      el(l)%id(k) = boundc(3)
      else if (coord(el(l)%n(k),2)==yi) then
      el(l)%id(k) = boundc(2)
      else if (coord(el(l)%n(k),2)==yl) then
      el(l)%id(k) = boundc(4)
      else
      el(l)%id(k) = 0
      end if
!

!     Impor que nos do canto da malha sejam sempre de Dirichlet nos encontros 2 e 1
!
!      if (coord(el(l)%n(k),1)==xi .and. coord(el(l)%n(k),2)==yi) then
!      if (boundc(1)==1 .or. boundc(2)==1) el(l)%id(k) = 1
!
!      else if(coord(el(l)%n(k),1)==xl .and. coord(el(l)%n(k),2)==yi) then
!      if (boundc(2)==1 .or. boundc(3)==1) el(l)%id(k) = 1
!
!      else if(coord(el(l)%n(k),1)==xl .and. coord(el(l)%n(k),2)==yl) then
!      if (boundc(3)==1 .or. boundc(4)==1) el(l)%id(k) = 1
!
!      else if(coord(el(l)%n(k),1)==xi .and. coord(el(l)%n(k),2)==yl) then
!      if (boundc(4)==1 .or. boundc(1)==1) el(l)%id(k) = 1
!
!      end if
!
     um4 = xi+(xl-xi)/4.d0
     um2 = yi+(yl-yi)/2.d0
!
!....Degrau
!
!     if(coord(el(l)%n(k),1)<=um4.and.coord(el(l)%n(k),2)<=um2) then
!     el(l)%id(k) = 1
!     endif
!
!.....Ventriculo
!
!     if(coord(el(l)%n(k),1)>(xi+3.d0*(xl-xi)/4.d0).and.coord(el(l)%n(k),1)<xl.and.coord(el(l)%n(k),2)==yi)  el(l)%id(k) = 5
!
!...........
!
!      print *,  el(l)%id(k)
      end do ! k
      end do ! i
!
      end do ! j
!
!      do i=1,nelem
!      do j=1, nen
!      write(*,*) i, j, el(i)%id(j)
!      end do !j
!      end do !i
!      stop
!
!...................................
!
      case default
!
      write(*,*) 'element degree not defined'
      stop
!
      end select
!
!     conectivities: element-edges
!
!     vertical edges
!
      k=0
      do j=1,ny
!
      do i=1,nx
!
      l=i+(nx)*(j-1)
      k=k+1
      el(l)%e(1)=k
      k=k+1
      el(l)%e(3)=k
      k=k-1
!
      end do ! i
!
      k=k+1
!
      end do ! j
!
!     horizontal edges
!
      do i=1,nx
!
      do j=1,ny
!
      l=i+(nx)*(j-1)
      k=k+1
      el(l)%e(2)=k
      k=k+1
      el(l)%e(4)=k
      k=k-1
!
      end do ! i
!
      k=k+1
!
      end do ! j
!
!     neighbors
!
      do i=1,nx
      do j=1,ny
!
      l =i+(nx)*(j-1)
      l1=l-1
      l2=l-nx
      l3=l+1
      l4=l+nx
!
      if(i.eq.1 ) l1=0
      if(j.eq.1 ) l2=0
      if(i.eq.nx) l3=0
      if(j.eq.ny) l4=0
!
      el(l)%ngb(1)=l1
      el(l)%ngb(2)=l2
      el(l)%ngb(3)=l3
      el(l)%ngb(4)=l4
!
      end do ! i
      end do ! j
!
      end subroutine
!
!---------------------------------------------------------------------------
!
      subroutine gausstable
!
!     store data for numerical integration
!
      use mgeometry
      use mgauss
!
      implicit none
!
      integer :: i,j,k,l,nint1d,nint2d
      integer, parameter :: ngt=8
!      real(8), dimension(ngt,ngt) :: xig,wgl
!
!     unidimensional
!
!     one point
!
      xig(1,1) = 0.d0
!
      wgl(1,1) = 2.d0
!
!     two points
!
      xig(1,2) = -dsqrt(3.d0)/3.d0
      xig(2,2) =  dsqrt(3.d0)/3.d0
!
      wgl(1,2) =  1.d0
      wgl(2,2) =  1.d0
!
!     three points
!
      xig(1,3) =  0.d0
      xig(2,3) = -dsqrt(3.d0/5.d0)
      xig(3,3) =  dsqrt(3.d0/5.d0)
!
      wgl(1,3) =  8.d0/9.d0
      wgl(2,3) =  5.d0/9.d0
      wgl(3,3) =  5.d0/9.d0
!
!     four points
!
      xig(1,4) = -dsqrt((3.d0-2.d0*dsqrt(6.d0/5.d0))/7.d0)
      xig(2,4) =  dsqrt((3.d0-2.d0*dsqrt(6.d0/5.d0))/7.d0)
      xig(3,4) = -dsqrt((3.d0+2.d0*dsqrt(6.d0/5.d0))/7.d0)
      xig(4,4) =  dsqrt((3.d0+2.d0*dsqrt(6.d0/5.d0))/7.d0)
!
      wgl(1,4) =  (18.d0+dsqrt(30.d0))/36.d0
      wgl(2,4) =  (18.d0+dsqrt(30.d0))/36.d0
      wgl(3,4) =  (18.d0-dsqrt(30.d0))/36.d0
      wgl(4,4) =  (18.d0-dsqrt(30.d0))/36.d0
!
!     five points
!
      xig(1,5) =  0.d0
      xig(2,5) = -dsqrt(5.d0-2.d0*dsqrt(10.d0/7.d0))/3.d0
      xig(3,5) =  dsqrt(5.d0-2.d0*dsqrt(10.d0/7.d0))/3.d0
      xig(4,5) = -dsqrt(5.d0+2.d0*dsqrt(10.d0/7.d0))/3.d0
      xig(5,5) =  dsqrt(5.d0+2.d0*dsqrt(10.d0/7.d0))/3.d0
!
      wgl(1,5) =  128.d0/225.d0
      wgl(2,5) =  (322.d0+13.d0*dsqrt(70.d0))/900.d0
      wgl(3,5) =  (322.d0+13.d0*dsqrt(70.d0))/900.d0
      wgl(4,5) =  (322.d0-13.d0*dsqrt(70.d0))/900.d0
      wgl(5,5) =  (322.d0-13.d0*dsqrt(70.d0))/900.d0
!
!     six points
!
      xig(1,6) =  0.6612093864662645
      xig(2,6) = -0.6612093864662645
      xig(3,6) = -0.2386191860831969
      xig(4,6) =  0.2386191860831969
      xig(5,6) = -0.9324695142031521
      xig(6,6) =  0.9324695142031521
!   
      wgl(1,6) =  0.3607615730481386
      wgl(2,6) =  0.3607615730481386
      wgl(3,6) =  0.4679139345726910
      wgl(4,6) =  0.4679139345726910
      wgl(5,6) =  0.1713244923791704
      wgl(6,6) =  0.1713244923791704	
!
!     seven points
!
      xig(1,7) =  0.0000000000000000
      xig(2,7) =  0.4058451513773972
      xig(3,7) = -0.4058451513773972
      xig(4,7) = -0.7415311855993945
      xig(5,7) =  0.7415311855993945
      xig(6,7) = -0.9491079123427585
      xig(7,7) =  0.9491079123427585
!   	
      wgl(1,7) =  0.4179591836734694	
      wgl(2,7) =  0.3818300505051189	
      wgl(3,7) =  0.3818300505051189	
      wgl(4,7) =  0.2797053914892766	
      wgl(5,7) =  0.2797053914892766	
      wgl(6,7) =  0.1294849661688697	
      wgl(7,7) =  0.1294849661688697	
!
!     eight points      
!

      xig(1,8) = -0.1834346424956498
      xig(2,8) =  0.1834346424956498
      xig(3,8) = -0.5255324099163290
      xig(4,8) =  0.5255324099163290
      xig(5,8) = -0.7966664774136267
      xig(6,8) =  0.7966664774136267
      xig(7,8) = -0.9602898564975363
      xig(8,8) =  0.9602898564975363
!   	
      wgl(1,8) =  0.3626837833783620	
      wgl(2,8) =  0.3626837833783620
      wgl(3,8) =  0.3137066458778873	
      wgl(4,8) =  0.3137066458778873	
      wgl(5,8) =  0.2223810344533745	
      wgl(6,8) =  0.2223810344533745	
      wgl(7,8) =  0.1012285362903763	
      wgl(8,8) =  0.1012285362903763
!
!     bidimensional
!
      do nint1d=1,ngt
!
      nint2d=nint1d**nsd
!
      l=0
!
      do j=1,nint1d
!
      do i=1,nint1d
!
      l=l+1
!
      xi(l,nint2d,1) = xig(i,nint1d)
      xi(l,nint2d,2) = xig(j,nint1d)
!
      wg(l,nint2d)=wgl(i,nint1d)*wgl(j,nint1d)
!
      end do ! i
!
      end do ! j
!
      end do ! nint1d
!
!      do i=1,4
!      write(*,*) (xi(i,4,j),j=1,nsd)
!      end do ! i
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine assemble
!
      use mgeometry
      use mcoeficientes
      use mgauss
      use marrays
      use mtime
!
      implicit none
!
      integer :: nel,i,j,k,l,m
      integer :: nint=64,nint1d=8
      real(8), dimension(nsd,nsd) :: jf,invjf,tinvjf     ! matriz jacobiana, matriz inversa do jacobiano, matriz inversa transposta
      real(8), dimension(nsd)     :: xil,xxl             !
      real(8), dimension(nen,nsd) :: xe                  ! coord de x na malha real
      real(8), dimension(2*nen,2*nen) :: ae              !matriz de rigidez local do problema de elasticidade
      real(8), dimension(2*nen, np) :: be                ! matriz local correspondente a termos de pressão e velocidade
      real(8), dimension(np, 2*nen) :: de                ! matriz local correspondente a termos de pressão e velocidade
      real(8), dimension(np, np) :: ce                   ! matriz local correspondente aos termos de pressao>invmb> sua inversa
      real(8), dimension(2*nen,  1) :: fe                ! vetor fonte local do problema de elasticidade
      real(8), dimension(np,  1) :: ge                   ! vetor fonte da pressao locais
      real(8) :: detjf                                   ! determinante da jacobiana
      real(8) :: phi2d,dphi2dx,dphi2dy,psig              ! funcoes base das variaveis no elemento de referencia, suas derivadas em x e y
      real(8) :: dpsigx,dpsigy
      real(8) :: thi                                     ! funcao base da pressao das variaveis no elemento de referencia
      real(8) :: f,g,h,a                                 ! funcoes das variaveis na malha real
      integer :: ig,jg,ibanda,jbanda,ib1, ib2, ib3, ib4,jb1, jb2, jb3, jb4
      integer :: i1,i2,j1,j2
      real(8) :: d1,d2,dxi,dxj,dyi,dyj
      real(8) :: dhx, dhy                                ! derivadas parciais 
      real(8) :: dhdx, dhdy, dht                         ! derivadas totais
      real(8) :: area,hh
      real(8), dimension(nsd) :: uk,ua
      real(8) :: dpx,dpy,dhp
      real(8) :: L0, ta,tnp
      real(8), dimension(2):: teste
      real(8), dimension(nsd,nsd):: duk
      real(8), dimension(nsd):: duu
      real(8) :: xil1d
      real(8) :: delta1
      real(8), dimension(4)           :: ee
!
      flowr=0.d0
      ta=tt-dt
      if(tt==ti) ta=ti
!
!....Calculating the volume change frame to frame
!
!     do nel=1,nelem !faz o loop em todos os elementos
!
!     coordinates of the element nodes
!
!      do j=1,nen
!      do k=1,nsd
!      xe(j,k) = coord(el(nel)%n(j),k)
!      end do ! k
!      end do ! j
!
!      do l=1,nint
!
!     integration points
!
!      do k=1,nsd
!      xil(k) = xi(l,nint,k)
!      end do ! k
!
!     real coordinates of the integration point
!
!      xxl=0.d0
!
!      do j=1,neng
!      do k=1,nsd
!      xxl(k) = xxl(k) + phi2d(j,xil,neng)*xe(j,k)
!      end do ! k
!      end do ! j
!
!      call jacobf(xe,xil,jf,detjf)
!
!     calcula p no ponto de integracao
!
!      pp =0.d0
!
!      do m=1,np
!      pp =pp +psig(m,xxl)*cpres(m,nel)
!      end do ! m
!
!     Calculating the volume change frame to frame
!
!      flowr = flowr + dht(xxl,pp,tt)*detjf*wg(l,nint) !(h(xxl,pp,tt)-h(xxl,ppa,ta))
!
!      end do !l
!      end do !nel
!............................
!
!     write(*,*) flowr
!
!
      tnp=tt
!      tnp=tt-dt/2.d0
!
      do nel=1,nelem !faz o loop em todos os elementos
!
      ae = 0.d0
      be = 0.d0
      ce = 0.d0
      de = 0.d0
      fe = 0.d0
      ge = 0.d0
!
!     coordinates of the element nodes
!
      do j=1,nen
      do k=1,nsd
      xe(j,k) = coord(el(nel)%n(j),k)
!      write(*,*) nel,xe(j,k)
      end do ! k
      end do ! j
!
      ee(1)=dsqrt((xe(1,1)-xe(4,1))**2.d0+(xe(1,2)-xe(4,2))**2.d0)
      ee(2)=dsqrt((xe(2,1)-xe(1,1))**2.d0+(xe(2,2)-xe(1,2))**2.d0)
      ee(3)=dsqrt((xe(3,1)-xe(2,1))**2.d0+(xe(3,2)-xe(2,2))**2.d0)
      ee(4)=dsqrt((xe(4,1)-xe(3,1))**2.d0+(xe(4,2)-xe(3,2))**2.d0)
!
!     compute the diameter of the element
!
      area=0.d0
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k)=xi(l,nint,k)
      end do ! k
!
      call jacobf(xe,xil,jf,detjf)
!
      area=area+detjf*wg(l,nint)
!
      end do ! l
!
      hh=dsqrt(2.d0*area)
!
!      estabilizacao
!
      delta1=0.d0
!      if(tt.ge.t4) delta1=(10.d0/mu)*hh
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k) = xi(l,nint,k)
      end do ! k
!
!     real coordinates of the integration point
!
      xxl=0.d0
!
      do j=1,neng
      do k=1,nsd
      xxl(k) = xxl(k) + phi2d(j,xil,neng)*xe(j,k)
      end do ! k
      end do ! j
!
      call jacobf(xe,xil,jf,detjf)
!
!     Inverse of the Jacobian
!
      invjf(1,1) = jf(2,2)
      invjf(1,2) =-jf(1,2)
      invjf(2,1) =-jf(2,1)
      invjf(2,2) = jf(1,1)
!
      invjf = invjf/detjf
!
!     Transpose
!
      tinvjf = transpose(invjf)
!
!     calculando a solucao no instante anterior no ponto de integracao
!
      ua=0.d0
!
      do i=1,nen
      do k=1,nsd
      ua(k)=ua(k)+phi2d(i,xil,nen)*pa(nel,(i-1)*nsd+k)
      end do ! k
      end do !i
!
!     calculando a solucao na iteracao k no ponto de integracao
!
      uk=0.d0
!
      do i=1,nen
      do k=1,nsd
      uk(k)=uk(k)+phi2d(i,xil,nen)*pk(nel,(i-1)*nsd+k)
      end do ! k
      end do !i
!
!     calculando o gradiente da solucao na iteracao k no ponto de integracao
!
      duk=0.d0
!
      do i=1,nen
!      
      i1 = (i-1)*nsd + 1
      i2 = (i-1)*nsd + 2
!
      d1 = dphi2dx(i,xil,nen)
      d2 = dphi2dy(i,xil,nen)
!
      dxi = tinvjf(1,1)*d1 + tinvjf(1,2)*d2
      dyi = tinvjf(2,1)*d1 + tinvjf(2,2)*d2
!      
      duk(1,1)=duk(1,1)+dxi*p(nel,(i-1)*nsd+1)
      duk(1,2)=duk(1,2)+dyi*p(nel,(i-1)*nsd+1)
      duk(2,1)=duk(2,1)+dxi*p(nel,(i-1)*nsd+2)
      duk(2,2)=duk(2,2)+dyi*p(nel,(i-1)*nsd+2)
!
      end do !i
!
!     grad de u aplicado a u no instante anterior
!
      duu(1)=duk(1,1)*uk(1)+duk(1,2)*uk(2)
      duu(2)=duk(2,1)*uk(1)+duk(2,2)*uk(2)
!
!     calcula p e grad de p em k avaliado no ponto de integracao
!
      pp =0.d0
      ppa=0.d0
      dpx=0.d0
      dpy=0.d0
      do m=1,np
      pp =pp +psig  (m,xxl)*cpres (m,nel)
      ppa=ppa+psig  (m,xxl)*cpresa(m,nel)
      dpx=dpx+dpsigx(m,xxl)*cpres (m,nel)
      dpy=dpy+dpsigy(m,xxl)*cpres (m,nel)
      end do ! m
!
!     RIGIDEZ LOCAL
!
      do i=1,nen
!
      i1 = (i-1)*nsd + 1
      i2 = (i-1)*nsd + 2
!
      d1 = dphi2dx(i,xil,nen)
      d2 = dphi2dy(i,xil,nen)
!
      dxi = tinvjf(1,1)*d1 + tinvjf(1,2)*d2
      dyi = tinvjf(2,1)*d1 + tinvjf(2,2)*d2
!
      do j=1,nen
!
      j1 = (j-1)*nsd+1   
      j2 = (j-1)*nsd+2
!
      d1 = dphi2dx(j,xil,nen)
      d2 = dphi2dy(j,xil,nen)
!
      dxj = tinvjf(1,1)*d1 + tinvjf(1,2)*d2
      dyj = tinvjf(2,1)*d1 + tinvjf(2,2)*d2
!
!      dhdx=dhx(xxl,pp,tt)+dhp(xxl)*dpx
!      dhdy=dhy(xxl,pp,tt)+dhp(xxl)*dpy
!
!     TERMOS CLASSICOS (independem de h)
!
!
!     termo de ordem zero (matriz de massa) ATENCAO: h nao ta variando com p. tem que calcular p medio
!
      ae(i1,j1) = ae(i1,j1) + (h(xxl,pp,tnp)*rho)*phi2d(i,xil,nen)*phi2d(j,xil,nen)*detjf*wg(l,nint)
      ae(i2,j2) = ae(i2,j2) + (h(xxl,pp,tnp)*rho)*phi2d(i,xil,nen)*phi2d(j,xil,nen)*detjf*wg(l,nint)
!
!     termo de primeira ordem (adveccao)
!
!      ae(i1,j1) = ae(i1,j1) + dt*h(xxl,pp,tt)*rho*phi2d(j,xil,nen)*(uk(1)*dxi+uk(2)*dyi)*detjf*wg(l,nint)
!      ae(i2,j2) = ae(i2,j2) + dt*h(xxl,pp,tt)*rho*phi2d(j,xil,nen)*(uk(1)*dxi+uk(2)*dyi)*detjf*wg(l,nint)
!
!     termo de primeira ordem (adveccao)
!
      ae(i1,j1) = ae(i1,j1) + dt*h(xxl,pp,tt)*rho*phi2d(i,xil,nen)*(uk(1)*dxj+uk(2)*dyj)*detjf*wg(l,nint)
      ae(i2,j2) = ae(i2,j2) + dt*h(xxl,pp,tt)*rho*phi2d(i,xil,nen)*(uk(1)*dxj+uk(2)*dyj)*detjf*wg(l,nint)
!
!     usando a solucao no instante anterior
!
!      ae(i1,j1) = ae(i1,j1) + dt*h(xxl,pp,tt-dt/2.d0)*rho*phi2d(j,xil,nen)*(ua(1)*dxi+ua(2)*dyi)*detjf*wg(l,nint)
!      ae(i2,j2) = ae(i2,j2) + dt*h(xxl,pp,tt-dt/2.d0)*rho*phi2d(j,xil,nen)*(ua(1)*dxi+ua(2)*dyi)*detjf*wg(l,nint)
!
!      usando o gradiente na iteracao anterior
!
!      ae(i1,j1) = ae(i1,j1) + dt*h(xxl,pp,tt)*rho*(duk(1,1))*phi2d(i,xil,nen)*phi2d(j,xil,nen)*detjf*wg(l,nint)
!      ae(i2,j1) = ae(i2,j1) + dt*h(xxl,pp,tt)*rho*(duk(1,2))*phi2d(i,xil,nen)*phi2d(j,xil,nen)*detjf*wg(l,nint)
!
!      ae(i1,j2) = ae(i1,j2) + dt*h(xxl,pp,tt)*rho*(duk(2,1))*phi2d(i,xil,nen)*phi2d(j,xil,nen)*detjf*wg(l,nint)
!      ae(i2,j2) = ae(i2,j2) + dt*h(xxl,pp,tt)*rho*(duk(2,2))*phi2d(i,xil,nen)*phi2d(j,xil,nen)*detjf*wg(l,nint)
!
!     Termo de segunda ordem (viscosidade)
!
!     primeira parte
!
!      ae(i1,j1) = ae(i1,j1) + dt*(mu*(1.d0+delta1)*(2.d0*dxi*dxj+dyi*dyj))*detjf*wg(l,nint)
!      ae(i2,j2) = ae(i2,j2) + dt*(mu*(1.d0+delta1)*(dxi*dxj+2.d0*dyi*dyj))*detjf*wg(l,nint)
!
!     segunda parte 
!
!      ae(i1,j2) = ae(i1,j2) + dt*mu*(1.d0+delta1)*(dyi*dxj)*detjf*wg(l,nint)
!      ae(i2,j1) = ae(i2,j1) + dt*mu*(1.d0+delta1)*(dxi*dyj)*detjf*wg(l,nint)
!
!     primeira parte
!
      ae(i1,j1) = ae(i1,j1) + dt*(mu*(1.d0+delta1)*(2.d0*dxi*dxj+dyi*dyj))*detjf*wg(l,nint)
      ae(i2,j2) = ae(i2,j2) + dt*(mu*(1.d0+delta1)*(dxi*dxj+2.d0*dyi*dyj))*detjf*wg(l,nint)
!
!     segunda parte 
!
      ae(i1,j2) = ae(i1,j2) + dt*mu*(1.d0+delta1)*(dyj*dxi)*detjf*wg(l,nint)
      ae(i2,j1) = ae(i2,j1) + dt*mu*(1.d0+delta1)*(dxj*dyi)*detjf*wg(l,nint)
!
      end do ! j
!
!     matriz B
!
      do m=1, np
!
!     TERMOS CLASSICOS (independem de h)
!
      be(i1,m) = be(i1,m) - dt*dxi*psig(m,xxl)*detjf*wg(l,nint)
      be(i2,m) = be(i2,m) - dt*dyi*psig(m,xxl)*detjf*wg(l,nint)
!
      end do !m
!
!     termo de fonte da equacao do equilibrio !coloquei em tn+1
!
      fe(i1,1) = fe(i1,1) + (dt*h(xxl,pp,tt)*rho)*phi2d(i,xil,nen)*f(xxl,1)*detjf*wg(l,nint)
      fe(i2,1) = fe(i2,1) + (dt*h(xxl,pp,tt)*rho)*phi2d(i,xil,nen)*f(xxl,2)*detjf*wg(l,nint)
!
!     adicionando os termos da solucao no instante anterior
!
      fe(i1,1) = fe(i1,1)+(h(xxl,pp,tnp)*rho)*ua(1)*phi2d(i,xil,nen)*detjf*wg(l,nint)
      fe(i2,1) = fe(i2,1)+(h(xxl,pp,tnp)*rho)*ua(2)*phi2d(i,xil,nen)*detjf*wg(l,nint)
!
!      jogando a advecao pro termo de fonte
! 
!      fe(i1,1) = fe(i1,1)-dt*h(xxl,pp,tt)*rho*phi2d(i,xil,nen)*duu(1)*detjf*wg(l,nint)
!      fe(i2,1) = fe(i2,1)-dt*h(xxl,pp,tt)*rho*phi2d(i,xil,nen)*duu(2)*detjf*wg(l,nint)
!
!     CONSERVACAO DE MASSA
!
      do m=1, np
!
!     matriz D !!!no tempo atual
!
      de(m,i1) = de(m,i1) + (h(xxl,pp,tt)*dxi+phi2d(i,xil,nen)*dhx(xxl,pp,tt))*psig(m,xxl)*detjf*wg(l,nint)
      de(m,i2) = de(m,i2) + (h(xxl,pp,tt)*dyi+phi2d(i,xil,nen)*dhy(xxl,pp,tt))*psig(m,xxl)*detjf*wg(l,nint)
!
      end do !m
!
      end do ! i
!
!     matriz C
!
      do m=1,np
!
      do j=1,np
!
      ce(m,j) = ce(m,j) + (1.d0/lambda)*psig(m,xxl)*psig(j,xxl)*detjf*wg(l,nint) !no limite nao importa o sinal
!
!      ce(m,j) = ce(m,j) - (1.d0/lambda)*psig(m,xil)*psig(j,xil)*detjf*wg(l,nint)
!
      end do !j
!
!     termo de fonte do balanco de massa
!
!      ge(m,1) = ge(m,1) + psig(m,xxl)*(g(xxl)*h(xxl,pp,tt)-(h(xxl,pp,tt)-h(xxl,ppa,ta))/dt)*detjf*wg(l,nint) !
      ge(m,1) = ge(m,1) + psig(m,xxl)*(g(xxl)*h(xxl,pp,tt)-dht(xxl,pp,tt))*detjf*wg(l,nint) !(h(xxl,pp,tt)-h(xxl,ppa,ta))/dt
!
!      if(l==nint)   write(*,*) ge(m,1)
!
!      ge(m,1) = ge(m,1) - psig(m,xil)*g(xxl)*h(xxl,pp,tt)*detjf*wg(l,nint)
!
      end do !m
!
      end do ! l
!
!     calcula as integrais na fronteira para fazer robin
!
      if(coord(el(nel)%n(7),2).eq.yint(2)) then
!
!      if(coord(el(nel)%n(7),1).ge.xai.and.coord(el(nel)%n(7),1).le.xaf.or. &
!         coord(el(nel)%n(7),1).ge.xmi.and.coord(el(nel)%n(7),1).le.xmf) then
!
      if(coord(el(nel)%n(7),1).ge.xai.and.coord(el(nel)%n(7),1).le.xaf) then
!
      do l=1,nint1d
!
      xil1d = xig(l,nint1d)
!
!     integration point
!
      xil(1) = xil1d
      xil(2) = 1.d0
!
!     real coordinates of the integration point
!
      xxl=0.d0
!
      do j=1,neng
      do k=1,nsd
      xxl(k) = xxl(k) + phi2d(j,xil,neng)*xe(j,k)
      end do ! k
      end do ! j
!
      call jacobf(xe,xil,jf,detjf)
!
!     Inverse of the Jacobian
!
      invjf(1,1) = jf(2,2)
      invjf(1,2) =-jf(1,2)
      invjf(2,1) =-jf(2,1)
      invjf(2,2) = jf(1,1)
!
      invjf = invjf/detjf
!
!     Transpose
!
      tinvjf = transpose(invjf)
!
      do i=1,nen
!
      if(i.eq.3.or.i.eq.4.or.i.eq.7) then
!
      i1 = (i-1)*nsd + 1
      i2 = (i-1)*nsd + 2
!
      d1 = dphi2dx(i,xil,nen)
      d2 = dphi2dy(i,xil,nen)
!
      dxi = tinvjf(1,1)*d1 + tinvjf(1,2)*d2
      dyi = tinvjf(2,1)*d1 + tinvjf(2,2)*d2
!
      do j=1,nen
!
      if(j.eq.3.or.j.eq.4.or.j.eq.7) then
!
      j1 = (j-1)*nsd+1   
      j2 = (j-1)*nsd+2
!
      d1 = dphi2dx(j,xil,nen)
      d2 = dphi2dy(j,xil,nen)
!
      dxj = tinvjf(1,1)*d1 + tinvjf(1,2)*d2
      dyj = tinvjf(2,1)*d1 + tinvjf(2,2)*d2
!
!     primeira parte
!
      ae(i1,j1) = ae(i1,j1) - dt*(     mu*dyj*phi2d(i,xil,nen))*ee(4)/2.d0*wgl(l,nint1d)
      ae(i2,j2) = ae(i2,j2) - dt*(2.d0*mu*dyj*phi2d(i,xil,nen))*ee(4)/2.d0*wgl(l,nint1d)
!
!     segunda parte 
!
      ae(i1,j2) = ae(i1,j2) - 0.d0
      ae(i2,j1) = ae(i2,j1) - dt*(     mu*dxj*phi2d(i,xil,nen))*ee(4)/2.d0*wgl(l,nint1d)
!
      end if
      end do ! j
!
!     matriz B
!
!      do m=1, np
!
!     parte da pressao
!
!      be(i1,m) = be(i1,m) + 0.d0
!      be(i2,m) = be(i2,m) + dt*phi2d(i,xil,nen)*psig(m,xxl)*ee(4)/2.d0*wgl(l,nint1d)
!
!      end do !m
!
      end if
      end do ! i
!
!      do m=1, np
!
!     matriz D !!!no tempo atual
!
!      de(m, 6) = de(m, 6) + dt*psig(m,xxl)*phi2d(3,xil,nen)*ee(4)/2.d0*wgl(l,nint1d)
!      de(m, 8) = de(m, 8) + dt*psig(m,xxl)*phi2d(4,xil,nen)*ee(4)/2.d0*wgl(l,nint1d)
!      de(m,14) = de(m,14) + dt*psig(m,xxl)*phi2d(7,xil,nen)*ee(4)/2.d0*wgl(l,nint1d)
!
!      end do !m
!
!      do m=1, np
!
!     matriz B
!
!      be( 6,m) = be( 6,m) + dt*psig(m,xxl)*phi2d(3,xil,nen)*ee(4)/2.d0*wgl(l,nint1d)
!      be( 8,m) = be( 8,m) + dt*psig(m,xxl)*phi2d(4,xil,nen)*ee(4)/2.d0*wgl(l,nint1d)
!      be(14,m) = be(14,m) + dt*psig(m,xxl)*phi2d(7,xil,nen)*ee(4)/2.d0*wgl(l,nint1d)
!
!      end do !m
!
      end do ! l
!
      end if
!
      end if
!
      call invmb(ce,np) !a partir daqui ce eh inversa de ce
!
!     ce = lambda*ce
!
!     save local matrices
!
      api(nel,:,:) = matmul(ce,de)
      ag (nel,:,:) = matmul(ce,ge)
!
!     montando sistema local ke*ue = he
!
      ae = ae - matmul(be,api(nel,:,:))
      fe = fe - matmul(be,ag(nel,:,:))
!
!      write(222,*) 'ae ',nel
!      do i=1,2*nen
!      write(222,'(100e15.8)') (ae(i,j),j=1,2*nen)
!      end do ! i
!      write(222,*)
!
      call contorno2(ae,de,fe,nel) 
!
!      write(222,*) 'ae depois ',nel
!      do i=1,2*nen
!      write(222,'(100e15.8)') (ae(i,j),j=1,2*nen)
!      end do ! i
!      write(222,*)
!
!     espalha a local na global
!
      do i=1,nen
!
      ig=el(nel)%n(i)
      bb(2*ig  ,1) = bb(2*ig  ,1) + fe(2*i  ,1)
      bb(2*ig-1,1) = bb(2*ig-1,1) + fe(2*i-1,1)
!
      do j=1,nen
!
      jg = el(nel)%n(j)
!
!     estrutura de banda para cada caso > lapack
!
      ib1 = ibanda(2*ig-1,2*jg-1)
      jb1 = jbanda(2*ig-1,2*jg-1)
      ab(ib1,jb1) = ab(ib1,jb1)+ae(2*i-1,2*j-1)
!
      ib2 = ibanda(2*ig,2*jg)
      jb2 = jbanda(2*ig,2*jg)
      ab(ib2,jb2) = ab(ib2,jb2)+ae(2*i,2*j)
!
      ib3 = ibanda(2*ig-1,2*jg)
      jb3 = jbanda(2*ig-1,2*jg) 
      ab(ib3,jb3) = ab(ib3,jb3)+ae(2*i-1,2*j)
!
      ib4 = ibanda(2*ig,2*jg-1)
      jb4 = jbanda(2*ig,2*jg-1)
      ab(ib4,jb4) = ab(ib4,jb4)+ae(2*i,2*j-1)
!
!
      end do ! j
!
      end do ! i
!
      end do ! nel
!
!     teste(1)=0.045d0
!     teste(2)=0.025d0
!
!      write(*,*) h(teste,0.d0,tt)
!
      end subroutine
!
!----------------------------------------------------------------------
!
!     caso vetorial
!
      subroutine contorno2(ae,de,fe,nel)
!
      use mgeometry
      use marrays
      use mgauss
      use mcoeficientes
      use mtime
!
      implicit none
!
      integer :: nint=64,nint1d=8 !antes: 25 e 5
      integer :: nel
      integer :: i, j, k, nk, l, n1, n2, i1, i2, j1, j2, nj, nj1, nj2
      integer, dimension(2,2) :: m
!      integer, dimension(1,1) :: lnint 
      real(8) :: xil1d, f2cont,h, dphi2dx, dphi2dy
      real(8) :: phi2d, d1, d2, dxi, dxj, dyi, dyj,dy
      real(8), dimension(2*nen,2*nen) :: ae, at 
      real(8), dimension(np, 2*nen)   :: de, dtr                
      real(8), dimension(2*nen,  1)   :: fe, ft 
      real(8), dimension(4)           :: ee
      real(8), dimension(nen,nsd)     :: xe
      real(8), dimension(nsd)         :: xil,xxl
      real(8) :: nor, pr
      real(8), dimension(nsd,nsd)     :: jf,invjf,tinvjf 
      real(8) :: detjf, psig  
      real(8) :: um4
      real(8) :: eps=1.d-5
      real(8) :: trac1,trac2
      real(8) :: L0
      real(8) :: pp0,pp1,pp2,pp3,dvdy,f1,f2,f3,hy
!
      um4 = xint(1) + (xint(2) - xint(1))/4.d0
!
      do i = 1, nen
!
      i1 = 0
      i2 = 0
      j1 = 0
      j2 = 0
      at = 0.d0
      ft = 0.d0
      dtr= 0.d0
      m=0
!
      i1 = (i-1)*nsd + 1
      i2 = (i-1)*nsd + 2
!
!     Impondo vertices da malha nulos
!
      if((nel==1.and.i==1).or.(nel==nx*ny.and.i==3).or.(nel==(nx*ny-nx+1).and.i==4).or.(nel==nx.and.i==2)) then
      fe(i1,1) = 0.d0
      fe(i2,1) = 0.d0
!
      do k = 1, 2*nen
!
      ae(i1,k ) = 0.d0
!      ae(k ,i1) = 0.d0
      ae(i2,k ) = 0.d0
!      ae(k ,i2) = 0.d0
!
      end do !k

      ae(i1,i1) = 1.d0
      ae(i2,i2) = 1.d0
      end if
!
!.....Integral de Tracao, t2=-p, t1=0
!
      if (el(nel)%id(i) == 2) then
!
!      print *, "entrou2", nel, i
!
!     Comprimento dos lados
!
      do j=1,nen
      do k=1,nsd
      xe(j,k) = coord(el(nel)%n(j),k)
      end do ! k
      end do ! j

      ee(1)=dsqrt((xe(1,1)-xe(4,1))**2.d0+(xe(1,2)-xe(4,2))**2.d0)
      ee(2)=dsqrt((xe(2,1)-xe(1,1))**2.d0+(xe(2,2)-xe(1,2))**2.d0)
      ee(3)=dsqrt((xe(3,1)-xe(2,1))**2.d0+(xe(3,2)-xe(2,2))**2.d0)
      ee(4)=dsqrt((xe(4,1)-xe(3,1))**2.d0+(xe(4,2)-xe(3,2))**2.d0)
!      
!      do nk = 1,nee
      if(coord(el(nel)%n(i),2).ge.(yint(2)-eps)) nk=4
      if(coord(el(nel)%n(i),2).le.(yint(1)+eps)) nk=2
      if(coord(el(nel)%n(i),1).ge.(xint(2)-eps)) nk=3
      if(coord(el(nel)%n(i),1).le.(xint(1)+eps)) nk=1
!
!     Calculo das integrais de linha
!
      do l=1,nint1d
!
      xil1d = xig(l,nint1d)
!
      xil(1) = xil1d
      xil(2) = xil1d
!
      if(nk.eq.1)      xil(1) = -1.d0
      if(nk.eq.2)      xil(2) = -1.d0
      if(nk.eq.3)      xil(1) =  1.d0
      if(nk.eq.4)      xil(2) =  1.d0
!
!     Real coordinates of the integration point
!
      xxl = 0.d0
!
      do n1=1,neng
      do n2=1,nsd
      xxl(n2) = xxl(n2) + phi2d(n1,xil,neng)*xe(n1,n2)
      end do ! n1
      end do ! n2
!
!     Vetor de tracao com pressao prescrita
!
      trac2=-pr(xxl)*nor(nk,2)
!
      fe(i2,1) = fe(i2,1) + dt*trac2*phi2d(i,xil,nen)*ee(nk)/2.d0*wgl(l,nint1d)
!      write(*,*) xxl(1),xxl(2),pr(xxl)
!
      end do ! l
!
 !     end do !nk
!
!      fe(i1,1) = f2cont(coord(el(nel)%n(i),1:nsd),1)
!
!      do k = 1, 2*nen
!
!      ae(i1,k ) = 0.d0
!      ae(k ,i1) = 0.d0
!
!      end do !k
!
!      ae(i1,i1) = 1.d0
!
!
!.....Integral de Tracao + Dirichlet // u1=0, t2=0
!
      elseif (el(nel)%id(i) == 3) then
!
!     Prescribing velocity in x
!
      fe(i1,1) = 0.d0
!      
      do k = 1, 2*nen
!
      ae(i1,k ) = 0.d0
!      ae(k ,i1) = 0.d0
!
      end do !k
!
      ae(i1,i1) = 1.d0
!
!     -------------------------
!
!.....Integral de Tracao + Dirichlet // u2=0, t1=0
!
      elseif (el(nel)%id(i) == 4) then
!
!     Prescribing velocity in y
!
      fe(i2,1) = f2cont(coord(el(nel)%n(i),1:nsd),2)
!      
      do k = 1, 2*nen
!
      ae(i2,k ) = 0.d0
!      ae(k ,i2) = 0.d0
!
      end do !k
!
      ae(i2,i2) = 1.d0
!
!     -------------------------

!.....Integral de Tracao + Dirichlet // u1=0, t2=-p
!
      elseif (el(nel)%id(i) == 5) then
!
!      print *, "entrou5", nel, i
!
!     Comprimento dos lados
!
      do j=1,nen
      do k=1,nsd
      xe(j,k) = coord(el(nel)%n(j),k)
      end do ! k
      end do ! j

      ee(1)=dsqrt((xe(1,1)-xe(4,1))**2.d0+(xe(1,2)-xe(4,2))**2.d0)
      ee(2)=dsqrt((xe(2,1)-xe(1,1))**2.d0+(xe(2,2)-xe(1,2))**2.d0)
      ee(3)=dsqrt((xe(3,1)-xe(2,1))**2.d0+(xe(3,2)-xe(2,2))**2.d0)
      ee(4)=dsqrt((xe(4,1)-xe(3,1))**2.d0+(xe(4,2)-xe(3,2))**2.d0)
!      
!      do nk = 1,nee
      if(coord(el(nel)%n(i),2).ge.(yint(2)-eps)) nk=4
      if(coord(el(nel)%n(i),2).le.(yint(1)+eps)) nk=2
      if(coord(el(nel)%n(i),1).ge.(xint(2)-eps)) nk=3
      if(coord(el(nel)%n(i),1).le.(xint(1)+eps)) nk=1
!
!     Calculo das integrais de linha
!
      do l=1,nint1d
!
      xil1d = xig(l,nint1d)
!
      xil(1) = xil1d
      xil(2) = xil1d
!
      if(nk.eq.1)      xil(1) = -1.d0
      if(nk.eq.2)      xil(2) = -1.d0
      if(nk.eq.3)      xil(1) =  1.d0
      if(nk.eq.4)      xil(2) =  1.d0
!
!     Real coordinates of the integration point
!
      xxl = 0.d0
!
      do n1=1,neng
      do n2=1,nsd
      xxl(n2) = xxl(n2) + phi2d(n1,xil,neng)*xe(n1,n2)
      end do ! n1
      end do ! n2
!
!     Vetor de tracao com pressao prescrita
!
!
      trac2=-pr(xxl)*nor(nk,2)
!
      fe(i2,1) = fe(i2,1) + dt*trac2*phi2d(i,xil,nen)*ee(nk)/2.d0*wgl(l,nint1d)
!      write(*,*) xxl(1),xxl(2),pr(xxl)
!
      end do ! l
!
 !     end do !nk

      fe(i1,1) = f2cont(coord(el(nel)%n(i),1:nsd),1)
!
      do k = 1, 2*nen
!
      ae(i1,k ) = 0.d0
!      ae(k ,i1) = 0.d0
!
      end do !k

      ae(i1,i1) = 1.d0
!
!.....Integral de Tracao + Dirichlet // u1=0, t2=-pressao iterativa
!
      elseif (el(nel)%id(i) == 6) then
!
!      print *, "entrou6",tt, nel, i
!
!     Comprimento dos lados
!
      do j=1,nen
      do k=1,nsd
      xe(j,k) = coord(el(nel)%n(j),k)
      end do ! k
      end do ! j

      ee(1)=dsqrt((xe(1,1)-xe(4,1))**2.d0+(xe(1,2)-xe(4,2))**2.d0)
      ee(2)=dsqrt((xe(2,1)-xe(1,1))**2.d0+(xe(2,2)-xe(1,2))**2.d0)
      ee(3)=dsqrt((xe(3,1)-xe(2,1))**2.d0+(xe(3,2)-xe(2,2))**2.d0)
      ee(4)=dsqrt((xe(4,1)-xe(3,1))**2.d0+(xe(4,2)-xe(3,2))**2.d0)
!      
!      do nk = 1,nee
      if(coord(el(nel)%n(i),2).ge.(yint(2)-eps)) nk=4
      if(coord(el(nel)%n(i),2).le.(yint(1)+eps)) nk=2
      if(coord(el(nel)%n(i),1).ge.(xint(2)-eps)) nk=3
      if(coord(el(nel)%n(i),1).le.(xint(1)+eps)) nk=1
!
!     Calculo das integrais de linha
!
      do l=1,nint1d
!
      xil1d = xig(l,nint1d)
!
      xil(1) = xil1d
      xil(2) = 1.d0
!
!     Real coordinates of the integration point
!
      xxl = 0.d0
!
      do n1=1,neng
      do n2=1,nsd
      xxl(n2) = xxl(n2) + phi2d(n1,xil,neng)*xe(n1,n2)
      end do ! n1
      end do ! n2
!
!     calcular a pressao no ponto de integracao
!
      pp1=pressao(nel,4)
      pp2=pressao(nel,7)
      pp3=pressao(nel,3)
      pp0=pp1*phi2d(4,xil,nen)+pp2*phi2d(7,xil,nen)+pp3*phi2d(3,xil,nen)
!      write(*,*) pp0,pr(xxl),pp2
!
!     fluxo normal no ponto de integracao
!
      hy=coord(el(nel)%n(4),2)-coord(el(nel)%n(1),2)
      f1=p(nel, 8)
      f2=p(nel,14)
      f3=p(nel, 6)
      dvdy=(f1*dphi2dy(4,xil,nen)+f2*dphi2dy(7,xil,nen)+f3*dphi2dy(3,xil,nen))*2.d0/hy
!
!     Vetor de tracao com pressao prescrita
!
!
!      trac2=0.d0!
      trac2=-pr(xxl)*nor(nk,2)
!      trac2=-pp0*nor(nk,2)+2.d0*mu*dvdy*nor(nk,2)
!      trac2=-pp0*nor(nk,2)!+2.d0*mu*dvdy*nor(nk,2)
!      trac2=-pr(xxl)*nor(nk,2)+2.d0*mu*dvdy*nor(nk,2)
!
      fe(i2,1) = fe(i2,1) + dt*trac2*phi2d(i,xil,nen)*ee(4)/2.d0*wgl(l,nint1d)
!      write(*,*) xxl(1),xxl(2),pr(xxl)
!
      end do ! l
!
 !     end do !nk
!
!      fe(i1,1) = 0.d0!f2cont(coord(el(nel)%n(i),1:nsd),1)
!
!      do k = 1, 2*nen
!
!      ae(i1,k ) = 0.d0
!      ae(k ,i1) = 0.d0
!
!      end do !k
!
!      ae(i1,i1) = 1.d0
!
      elseif (el(nel)%id(i) == 7) then
!
!
!     -------------------------
!
!.....Dirichlet
!    
      else if (el(nel)%id(i)==1) then
!
!      print *, "entrou1", nel, i
!
!      do j = 1, nen
!
!      j1 = (j-1)*nsd+1
!      j2 = (j-1)*nsd+2
!
!      if (el(nel)%id(j).ne.1) then
!
!      fe(j1,1) = fe(j1,1) - f2cont(coord(el(nel)%n(i),1:nsd),1)*ae(j1,i1) &
!                          - f2cont(coord(el(nel)%n(i),1:nsd),2)*ae(j1,i2)
!      fe(j2,1) = fe(j2,1) - f2cont(coord(el(nel)%n(i),1:nsd),1)*ae(j2,i1) &
!                          - f2cont(coord(el(nel)%n(i),1:nsd),2)*ae(j2,i2)
!      end if
!
!      end do !j
!
      fe(i1,1) = f2cont(coord(el(nel)%n(i),1:nsd),1)
      fe(i2,1) = f2cont(coord(el(nel)%n(i),1:nsd),2)
!
      do k = 1, 2*nen
!
      ae(i1,k ) = 0.d0
!      ae(k ,i1) = 0.d0
      ae(i2,k ) = 0.d0
!      ae(k ,i2) = 0.d0
!
      end do !k

      ae(i1,i1) = 1.d0
      ae(i2,i2) = 1.d0
!     -------------------------
!
      elseif(el(nel)%id(i)==0) then
!
!     inside mesh node
!
      else
!
      write(*,*) "Boundary type not defined", el(nel)%id(i)
!
      end if
!
!      if(nel==4)      write(59,*)  nel,i, fe(2*i-1,1),fe(2*i  ,1)
!
      end do !i
!
      end subroutine
!
!---------------------------------------------------------------
! 
!     coordenadas do vetor de fonte
!
      function f(xx,i)
      use mgeometry, only : nsd, bc,xint,yint
      use mcoeficientes
      use mtime, only: tt,ti,tf
      use marrays, only: pp
      implicit none
!
      real(8), dimension(nsd) :: xx
      real(8) :: hel
      real(8) :: f,pi,dpi,pi2, pe, dhx, dhy, h, pr
      integer :: i 
!
      pi =4.d0*datan(1.d0)
      dpi=2.d0*pi
      pi2 = pi*pi
!
      select case(i)
!
!.........
!
      case(1) 
!
      if(bc(1) == 1)      f = 1.d0/h(xx,pp,tt)*(-2.d0*mu*(xx(2)*dhy(xx,pp,tt)+h(xx,pp,tt)) &
                            +h(xx,pp,tt)*pi*dcos(pi*xx(1))*(xx(2)+1)*(xx(2)-1)+pr(xx)*dhx(xx,pp,tt))
!      if(bc(2) == 1)      f = 0.d0!1.d0/h(xx,pp,tt)*(xx(2)*dhy(xx,pp,tt)+pr(xx)*dhx(xx,pp,tt))
      if(bc(2) == 1)      f = 1.d0/h(xx,pp,tt)*((xx(2)-1.d0/2.d0)*dhy(xx,pp,tt)+pr(xx)*dhx(xx,pp,tt))
      if(bc(3) == 1)      f = pe(xx,1)*(mu*pi2-mu/4.d0/pi2-2.d0*pi2+(-2.d0*pi2*mu/h(xx,pp,tt)-mu/2.d0/h(xx,pp,tt))*dhy(xx,pp,tt))+&
			    & pe(xx,2)/h(xx,pp,tt)*(mu-1.d0)*dhx(xx,pp,tt)
      if(bc(4) == 1)      f = (2.d0*mu*pi2+1.d0+2.d0*pi2)/h(xx,pp,tt)*&
                            &(h(xx,pp,tt)*pe(xx,1)+dhx(xx,pp,tt)*dsin(pi*xx(1))*dsin(pi*xx(2)))
      if(bc(5) == 1)      f = 0.d0
      if(bc(6) == 1)then
      f=0.d0
      if(tt>0.05.and.tt<0.3) then
           f=0.d0!1.d1*(tt-ti)/(0.05-ti)*dcos(pi*((xx(1)-xint(1))/(xint(2)-xint(1)))) &
              !                     *dsin(pi*((xx(2)-yint(1))/(yint(2)-yint(1))))
      endif
      end if
!
!...........
!
       case(2)
!
       if(bc(1) == 1)      f = 1.d0/h(xx,pp,tt)*(-2.d0*mu*xx(2)*dhx(xx,pp,tt)+&
                             & h(xx,pp,tt)*2.d0*xx(2)*dsin(pi*xx(1))+pr(xx)*dhy(xx,pp,tt))
!       if(bc(2) == 1)      f = 0.d0!1.d0/h(xx,pp,tt)*(xx(2)*dhx(xx,pp,tt)+pr(xx)*dhy(xx,pp,tt))
       if(bc(2) == 1)      f = 1.d0/h(xx,pp,tt)*((xx(2)-1.d0/2.d0)*dhx(xx,pp,tt)+pr(xx)*dhy(xx,pp,tt))
       if(bc(3) == 1)      f = pe(xx,1)*(-2*pi2*mu)/h(xx,pp,tt)*(1.d0+1.d0/4.d0/pi2)*dhx(xx,pp,tt)+&
			     & pe(xx,2)*((-mu/4.d0-1.d0)*dhy(xx,pp,tt)/h(xx,pp,tt)-1.d0/2.d0+mu*pi2-mu/4.d0)
       if(bc(4) == 1)      f = (-2.d0*mu*pi2+1.d0+2.d0*pi2)/h(xx,pp,tt)&
                             &*(-h(xx,pp,tt)*pe(xx,2)+dhy(xx,pp,tt)*dsin(pi*xx(1))*dsin(pi*xx(2)))
       if(bc(5) == 1)      f = 0.d0
       if(bc(6) == 1) then
       f=0.d0
       if(tt>0.05.and.tt<0.3) then
       f=0.d0!1.d1*(tt-ti)/(0.05-ti)*dcos(pi*((xx(2)-yint(1))/(yint(2)-yint(1))))&
             !                  *dsin(pi*((xx(1)-xint(1))/(xint(2)-xint(1))))
       endif
       end if
!
       end select 
!
       end function
!
!----------------------------------------------------------------------
!
!     fonte do balanco de massa
!
      function g(xx)
      use mgeometry, only : nsd, bc
      use mcoeficientes
      implicit none
!
      real(8), dimension(nsd) :: xx
      real(8) :: g,pi,dpi,pi2, pr, h, dhy, dhx, pe, dpex, dpey
      integer :: i 
!
      pi =4.d0*datan(1.d0)
      dpi=2.d0*pi
      pi2 = pi*pi
!
!      if (bc(1) == 1)     g = 0.d0      
!      if (bc(2) == 1)     g = 1.d0/h(xx,pp,tt)*(pe(xx,1)*dhx(xx,pp,tt)+pe(xx,2)*dhy(xx,pp,tt))
!      if (bc(3) == 1)     g = pe(xx,1)*dhx(xx,pp,tt)+pe(xx,2)*dhy(xx,pp,tt)
!      if (bc(4) == 1)     g = pi/h(xx,pp,tt)*(dcos(pi*xx(1))*dsin(pi*xx(2))*dhx(xx,pp,tt)&
!
			   g = 0.d0
!                          g = dpex(xx,1)+dpey(xx,2) + 1.d0/h(xx,pp,tt)*(pe(xx,1)*dhx(xx,pp,tt)+pe(xx,2)*dhy(xx,pp,tt))
!
      end function
!
!----------------------------------------------------------------------
!
!     funcao espessura baseada na elastancia
!
      function hel(xx,pp)
      use mgeometry, only : nsd
      use mcoeficientes
      use mtime, only : tt
      implicit none
!
      real(8), dimension(nsd) :: xx
      real(8) :: pi,fluxtime,pp
      real(8) :: hel
      real(8) :: Emax, Emin, a, b, k0, k1, vd, th
      real(8) :: tce
      real(8) :: phi, E
!
      pi = 4.d0*datan(1.d0)
      Emax = 331972744.66d0
      Emin = 6532796.98d0
      a = 0.9d0
      b = 0.25d0
      k0 = 0.29d0
      k1 = 0.2d0
      Vd = 1.d-5
      th = 1.d0
!
      tce = k0+k1*th
!
      if(tt<=tce) phi = a*dsin(pi*tt/tce)-b*dsin(2*pi*tt/tce)
      if(tt>tce.and.tt<th) phi=0.d0
!
      E = Emin*(1-phi) + Emax*phi
!
      hel = (Vd-9.d-5 + pp/E)/14729.9d-5
!
      end function
!
!----------------------------------------------------------------------
!
!     funcao espessura
!
      function h(xx,pp,t)
      use mgeometry, only : nsd, xint, yint,h0
      use mcoeficientes
      implicit none
!
      real(8), dimension(nsd) :: xx, v4v, v6v
      real(8) :: pi,fluxtime,pp
      real(8) :: h,t, ly, lx,tau
      real(8) :: a,b,aex1,aex2,bex1,bex2
      real(8) :: t13,t1d13,t3d13, interpol1
!
      pi = 4.d0*datan(1.d0)
!
!     Antigo para malha 0.05 x 0.05
!
      if(t>0.05.and.t<=0.3) then
      h = 0.036d0+(2477.d0*xx(2)-155540.d0*xx(2)**2+9.d0*471111.d0*xx(2)**3 -(42400000.d0*xx(2)**4))/9.d0*&
         & (-0.09756d0*t+0.017073d0)
      elseif(t<=0.05) then
      h = 0.036d0+(2477.d0*xx(2)-155540.d0*xx(2)**2+9.d0*471111.d0*xx(2)**3 -(42400000.d0*xx(2)**4))/9.d0*&
         & (-0.09756d0*0.05+0.017073d0)
      elseif(t>0.3) then
      h = 0.036d0+(2477.d0*xx(2)-155540.d0*xx(2)**2+9.d0*471111.d0*xx(2)**3 -(42400000.d0*xx(2)**4))/9.d0*&
         & (-0.09756d0*0.3+0.017073d0)
      endif
!
!     Para qlq malha
!
      ly=yint(2)-yint(1)
      lx=xint(2)-xint(1)
!
      a=  4.d0*(vmin-vmax)/(lx*ly*1.20889)
      b=-0.7d0*(vmin-vmax)/(lx*ly*1.20889)
!
      t13=t3-t1
      t1d13=t1/t13
      t3d13=t3/t13
      v4v(1)=t4
      v4v(2)=v4
      v6v(1)=t6
      v6v(2)=v6
!      write(*,*) t13, t1d13, t3d13
!
      aex1= v3**(-t1d13)*v1**(t3d13)
      bex1= t13/dlog(v1/v3)
      aex2= vmax**(-0.45d0/0.55d0)*vmin**(1/0.55d0)
      bex2= (-0.55d0)/dlog(vmin/vmax)
!
      if(t<=t1)                tau=t1
      if(t>t1.and.t<=t3)       tau=t
      if(t>t3.and.t<=t4)       tau=t3
      if(t>t4.and.t<=t6)       tau=t
      if(t>t6) then
      t=t6
      write(*,*) "tempo maior que um ciclo 1", t
      endif
!
      if(t<=0.45) then
      h=h0+(32.d0*xx(2)/(3.d0*ly)-416.d0*xx(2)**2/(15.d0*ly**2)+512.d0*xx(2)**3/(15*ly**3)-256.d0*xx(2)**4/(15.d0*ly**4))*(a*tau+b)
      h=h0+(32.d0*xx(2)/(3.d0*ly)-416.d0*xx(2)**2/(15.d0*ly**2)+512.d0*xx(2)**3/(15*ly**3)-256.d0*xx(2)**4/(15.d0*ly**4))*&
       &(aex1*dexp(-tau/bex1)-vmed)/(1.20889d0*lx*ly)
      else
!      h=h0+(32.d0*xx(2)/(3.d0*ly)-416.d0*xx(2)**2/(15.d0*ly**2)+512.d0*xx(2)**3/(15*ly**3)-256.d0*xx(2)**4/(15.d0*ly**4))*&
!       &(aex2*dexp(tau/bex2)-vmed)/(1.20889d0*lx*ly)
      h=h0+(32.d0*xx(2)/(3.d0*ly)-416.d0*xx(2)**2/(15.d0*ly**2)+512.d0*xx(2)**3/(15*ly**3)-256.d0*xx(2)**4/(15.d0*ly**4))*&
       &(interpol1(t,v4v,v6v)-vmed)/(1.20889d0*lx*ly)
!
      endif
!
      end function
!
!----------------------------------------------------------------------
!
!     derivada de h em p
!
      function dhp(xx) !dhp(xx,pp)
      use mgeometry, only : nsd
      use mcoeficientes
      use mtime, only : tt
      implicit none
!
      real(8), dimension(nsd) :: xx
      real(8) :: dhp, pi,fluxtime,pp
!
      pi = 4.d0*datan(1.d0)
!
!      dhp = (1.d0+dsin(pi*xx(1))*dsin(pi*xx(2)))*betapres
!      dhp = dsin(pi*xx(1))*dsin(pi*xx(2))*betapres
!      dhp = betapres*(dsin(pi*xx(1))*dsin(pi*xx(2)))
      dhp = 0.d0!betapres*(xx(1)-xx(1)**2)*(xx(2)-xx(2)**2)
!
      end function
!
!----------------------------------------------------------------------
!
!     derivada de h em t
!
      function dht(xx, pp,t) 
      use mgeometry, only : nsd, xint, yint,h0
      use mcoeficientes
      implicit none
!
      real(8), dimension(nsd) :: xx, v4v, v6v
      real(8) :: pi,fluxtime,pp
      real(8) :: dht,t, ly, lx,tau
      real(8) :: a,b,aex1,aex2,bex1,bex2
      real(8) :: t13,t1d13,t3d13,interpol1dt 
!
      pi = 4.d0*datan(1.d0)
!
!     Antigo para malha 0.05 x 0.05
!
      if(t>0.05.and.t<=0.3) then
      dht = (2477.d0*xx(2)-155540.d0*xx(2)**2+9.d0*471111.d0*xx(2)**3 -(42400000.d0*xx(2)**4))/9.d0*(-0.09756d0)
      elseif(t<=0.05) then
      dht = 0.d0
      elseif(t>0.3) then
      dht = 0.0d0
      endif
!
!     Para qlq malha
!
      ly=yint(2)-yint(1)
      lx=xint(2)-xint(1)
!
      a=  4.d0*(vmin-vmax)/(lx*ly*1.20889)
      b=-0.7d0*(vmin-vmax)/(lx*ly*1.20889)
!
      t13=t3-t1
      t1d13=t1/t13
      t3d13=t3/t13
      v4v(1)=t4
      v4v(2)=v4
      v6v(1)=t6
      v6v(2)=v6
!
      aex1= v3**(-t1d13)*v1**(t3d13)
      bex1= t13/dlog(v1/v3)
      aex2= vmax**(-0.45d0/0.55d0)*vmin**(1/0.55d0)
      bex2= (-0.55d0)/dlog(vmin/vmax)
!
      if(t<=0.05)                  tau=0.05
      if(t>0.05.and.t<=0.30)       tau=t
      if(t>0.30.and.t<=0.45)       tau=0.3
      if(t>0.45.and.t<=1.00)       tau=t
      if(t>tcicle) then
      t=tcicle
      write(*,*) "tempo maior que um ciclo 2", t
      endif
!
      if(t<=0.45) then
      dht=(32.d0*xx(2)/(3.d0*ly)-416.d0*xx(2)**2/(15.d0*ly**2)+512.d0*xx(2)**3/(15*ly**3)-256.d0*xx(2)**4/(15.d0*ly**4))*(a)
      dht=(32.d0*xx(2)/(3.d0*ly)-416.d0*xx(2)**2/(15.d0*ly**2)+512.d0*xx(2)**3/(15*ly**3)-256.d0*xx(2)**4/(15.d0*ly**4))*&
       &(-aex1/bex1*dexp(-tau/bex1))/(1.20889d0*lx*ly)
      else
!      dht=(32.d0*xx(2)/(3.d0*ly)-416.d0*xx(2)**2/(15.d0*ly**2)+512.d0*xx(2)**3/(15*ly**3)-256.d0*xx(2)**4/(15.d0*ly**4))*&
!       &(aex2/bex2*dexp(tau/bex2))/(1.20889d0*lx*ly)
      dht=(32.d0*xx(2)/(3.d0*ly)-416.d0*xx(2)**2/(15.d0*ly**2)+512.d0*xx(2)**3/(15*ly**3)-256.d0*xx(2)**4/(15.d0*ly**4))*&
       &(interpol1dt(v4v,v6v))/(1.20889d0*lx*ly)
      endif
!
      end function
!
!----------------------------------------------------------------------
!
!     derivada de h em x
!
      function dhx(xx,pp,t) !dhx(xx,pp,tt)
      use mgeometry, only : nsd, xint, yint,h0
      use mcoeficientes
      implicit none
!
      real(8), dimension(nsd) :: xx
      real(8) :: pi,fluxtime,pp
      real(8) :: dhx,t, ly, lx,a,b,tau
!
      pi = 4.d0*datan(1.d0)
!
      dhx = 0.d0
!
      end function
!
!----------------------------------------------------------------------
!
!     derivada de h em y
!
      function dhy(xx,pp,t) !dhy(xx,pp,tt)
      use mgeometry, only : nsd, xint, yint,h0
      use mcoeficientes
      implicit none
!
      real(8), dimension(nsd) :: xx, v4v,v6v
      real(8) :: pi,fluxtime,pp
      real(8) :: dhy,t, ly, lx,tau
      real(8) :: a,b,aex1,aex2,bex1,bex2
      real(8) :: t13, t1d13, t3d13, interpol1
!
      pi = 4.d0*datan(1.d0)
!
!     Antigo para malha 0.05 x 0.05
!
      if(t>0.05.and.t<=0.3) then
      dhy = (2477.d0-2.d0*155540.d0*xx(2)+3.d0*9.d0*471111.d0*xx(2)**2 -(4.d0*42400000.d0*xx(2)**3))/9.d0*&
         & (-0.09756d0*t+0.017073d0)
      elseif(t<=0.05) then
      dhy = (2477.d0-2.d0*155540.d0*xx(2)+3.d0*9.d0*471111.d0*xx(2)**2 -(4.d0*42400000.d0*xx(2)**3))/9.d0*&
         & (-0.09756d0*0.05d0+0.017073d0)
      elseif(t>0.3) then
      dhy = (2477.d0-2.d0*155540.d0*xx(2)+3.d0*9.d0*471111.d0*xx(2)**2 -(4.d0*42400000.d0*xx(2)**3))/9.d0*&
         & (-0.09756d0*0.3d0+0.017073d0)
      endif
!
!     Para qlq malha
!
      ly=yint(2)-yint(1)
      lx=xint(2)-xint(1)
!
      a=  4.d0*(vmin-vmax)/(lx*ly*1.20889)
      b=-0.7d0*(vmin-vmax)/(lx*ly*1.20889)
!
      t13=t3-t1
      t1d13=t1/t13
      t3d13=t3/t13
      v4v(1)=t4
      v4v(2)=v4
      v6v(1)=t6
      v6v(2)=v6
!
      aex1= v3**(-t1d13)*v1**(t3d13)
      bex1= t13/dlog(v1/v3)
      aex2= vmax**(-0.45d0/0.55d0)*vmin**(1/0.55d0)
      bex2= (-0.55d0)/dlog(vmin/vmax)
!
      if(t<=0.05)                  tau=0.05
      if(t>0.05.and.t<=0.30)       tau=t
      if(t>0.30.and.t<=0.45)       tau=0.3
      if(t>0.45.and.t<=1.00)       tau=t
      if(t>tcicle) then
      t=tcicle
      write(*,*) "tempo maior que um ciclo 3", t
      endif
!
      if(t<=0.45) then
      dhy=(32.d0/(3.d0*ly)-(832.d0*xx(2))/(15.d0*ly**2)+(512.d0*xx(2)**2)/(5.d0*ly**3)-(1024.d0*xx(2)**3)/(15.d0*ly**4))*(a*tau+b)
      dhy=(32.d0/(3.d0*ly)-(832.d0*xx(2))/(15.d0*ly**2)+(512.d0*xx(2)**2)/(5.d0*ly**3)-(1024.d0*xx(2)**3)/(15.d0*ly**4))*&
          &(aex1*dexp(-tau/bex1)-vmed)/(1.20889d0*lx*ly)
      else
!      dhy=(32.d0/(3.d0*ly)-(832.d0*xx(2))/(15.d0*ly**2)+(512.d0*xx(2)**2)/(5.d0*ly**3)-(1024.d0*xx(2)**3)/(15.d0*ly**4))*&
!       &(aex2*dexp(tau/bex2)-vmed)/(1.20889d0*lx*ly)
      dhy=(32.d0/(3.d0*ly)-(832.d0*xx(2))/(15.d0*ly**2)+(512.d0*xx(2)**2)/(5.d0*ly**3)-(1024.d0*xx(2)**3)/(15.d0*ly**4))*&
         &(interpol1(t,v4v,v6v)-vmed)/(1.20889d0*lx*ly)
!
      endif
!
      end function
!
!----------------------------------------------------------------------
!
!     caso vetorial: pegar solucao exata u
!
      function f2cont(xx,i)
      use mgeometry, only: nsd, bc
      use mcoeficientes
      implicit none
!
      real(8), dimension(nsd) :: xx
      real(8) :: f2cont,pi,dpi,pi2, pe
      integer :: i
!
      pi  = 4.d0*datan(1.d0)
      dpi = 2.d0*pi
      pi2 = pi*pi
!
      select case(i)
!
      case(1)
      f2cont = pe(xx,1)
!
      case(2)
      f2cont = pe(xx,2)
!
      end select
!
      end function
!
!-----------------------------------------------------------------------
!
!     vetor normal
!
      function nor(nk,i)
      use mgeometry, only : nsd
      use mcoeficientes
      implicit none
!
      integer :: nk,i
      real(8) :: nor !coord normal vector
!
      select case(nk)
!
      case(1) ! (-1,0)
!
      select case(i)
      case(1)
      nor = -1.d0
      case(2)
      nor = 0.d0
      end select
!
!.......
!
      case(2) ! (0,-1)
!
      select case(i)
      case(1)
      nor = 0.d0
      case(2)
      nor = -1.d0
      end select
!
!......
!
      case(3) ! (1,0)
!
      select case(i)
      case(1)
      nor = 1.d0
      case(2)
      nor = 0.d0
      end select
!
!......
!
      case(4) ! (0,1)
!
      select case(i)
      case(1)
      nor = 0.d0
      case(2)
      nor = 1.d0 
      end select
!
!......
!
      case default
!
      write(*,*) 'degree not defined'
      stop
!
      end select
!
      end function
!
!----------------------------------------------------------------------
!
      function phi(i,nen,xi)
!
!     one dimensional lagrangian functions
!
      implicit none
!
      real(8) :: phi,xi
      integer :: i,nen
!
      select case(nen)
!
      case(1)
!
!     constant
!
      phi=1.d0
!
      case(2)
!
!     linear
!
      if(i.eq.1) then
      phi=0.5d0*(1.d0-xi)
      else
      if(i.eq.2) then
      phi=0.5d0*(1.d0+xi)
      else
      write(*,*) 'inexistent node'
      stop
      end if
      end if
!
      case(3)
!
!     quadratic
!
      if(i.eq.1) then
      phi= 0.5d0*xi*(xi-1.d0)
      else
      if(i.eq.2) then
      phi=(1.d0-xi)*(1.d0+xi)
      else
      if(i.eq.3) then
      phi= 0.5d0*xi*(xi+1.d0)
      else
      write(*,*) 'inexistent node'
      stop
      end if
      end if
      end if
!
      case(4)
!
!     cubic
!
      if(i.eq.1) then
      phi= (xi+1.d0/3.d0)*(xi-1.d0/3.d0)*(xi-1.d0)*(-9.d0/16.d0)
      else
      if(i.eq.2) then
      phi= (xi+1.d0)*(xi-1.d0/3.d0)*(xi-1.d0)*(27.d0/16.d0)
      else
      if(i.eq.3) then
      phi= (xi+1.d0)*(xi+1.d0/3.d0)*(xi-1.d0)*(-27.d0/16.d0)
      else
      if(i.eq.4) then
      phi= (xi+1.d0)*(xi+1.d0/3.d0)*(xi-1.d0/3.d0)*(9.d0/16.d0)
      else
      write(*,*) 'inexistent node'
      stop
      end if
      end if
      end if
      end if
!
      case default
!
      write(*,*) 'degree not defined'
      stop
!
      end select
!
      end function phi
!
!-----------------------------------------------------------------------
!
      function dphi(i,nen,xi)
!
!     derivative of the one dimensional lagrangian functions
!
      implicit none
!
      real(8) :: dphi,xi
      integer :: i,nen
!
      select case(nen)
!
      case(1)
!
!     constant
!
      dphi=0.d0
!
      case(2)
!
!     linear
!
      if(i.eq.1) then
      dphi=-0.5d0
      else
      dphi= 0.5d0
      end if
!
      case(3)
!
!     quadratic
!
      if(i.eq.1) then
      dphi= xi-0.5d0
      else
      if(i.eq.2) then
      dphi=-2.d0*xi
      else
      dphi= xi+0.5d0
      end if
      end if
!
      case(4)
!
!     cubic
!
      if(i.eq.1) then
      dphi= (-27.d0*xi**2+18.d0*xi+1.d0)/16.d0
      else
      if(i.eq.2) then
      dphi= (3.d0*xi**2-2.d0/3.d0*xi-1.d0)*(27.d0/16.d0)
      else
      if(i.eq.3) then
      dphi= (3.d0*xi**2+2.d0/3.d0*xi-1.d0)*(-27.d0/16.d0)
      else
      dphi= (27.d0*xi**2+18.d0*xi-1.d0)/16.d0
      end if
      end if
      end if
!
      case default
!
      write(*,*) 'degree not defined'
      stop
!
      end select
!
      end function
!
!-----------------------------------------------------------------------
!
      function phi2d(i,xi,n)
      use mgeometry, only : nsd
      implicit none
!
      real(8) :: phi2d
      real(8) :: phi
      real(8), dimension(nsd) :: xi
      integer :: i,n1d,n
!
      select case(n)
!
!........................
!
      case(1)
!
!     constant
!
      phi2d=1.d0
!
!........................
!
      case(4)
!
!     bilinear
!
      n1d=2
!
      select case(i)
!
      case(1)
      phi2d=phi(1,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(2)
      phi2d=phi(2,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(3)
      phi2d=phi(2,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(4)
      phi2d=phi(1,n1d,xi(1))*phi(2,n1d,xi(2))
!

      case default
!
      write(*,*) '1 node not defined'
      stop
!
      end select
!
!........................
!
      case(9)
!
!     biquadratic
!
      n1d=3
!
      select case(i)
!
      case(1)
      phi2d=phi(1,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(2)
      phi2d=phi(3,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(3)
      phi2d=phi(3,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(4)
      phi2d=phi(1,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(5)
      phi2d=phi(2,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(6)
      phi2d=phi(3,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(7)
      phi2d=phi(2,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(8)
      phi2d=phi(1,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(9)
      phi2d=phi(2,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case default
!
      write(*,*) '2 node not defined'
      stop
!
      end select
!
!........................
!
      case(16)
!
!     bicubic
!
      n1d=4
!
      select case(i)
!
      case(1)
      phi2d=phi(1,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(2)
      phi2d=phi(4,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(3)
      phi2d=phi(4,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(4)
      phi2d=phi(1,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(5)
      phi2d=phi(2,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(6)
      phi2d=phi(3,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(7)
      phi2d=phi(4,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(8)
      phi2d=phi(4,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(9)
      phi2d=phi(3,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(10)
      phi2d=phi(2,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(11)
      phi2d=phi(1,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(12)
      phi2d=phi(1,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(13)
      phi2d=phi(2,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(14)
      phi2d=phi(3,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(15)
      phi2d=phi(3,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(16)
      phi2d=phi(2,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case default
!
      write(*,*) '3 node not defined'
      stop
!
      end select
!
      case default
!
      write(*,*) 'degree not defined'
      stop
!
      end select
!
      end function
!
!-----------------------------------------------------------------------
!
      function dphi2dx(i,xi,n)
      use mgeometry, only : nsd
      implicit none
!
      real(8) :: dphi2dx
      real(8) :: phi,dphi
      real(8), dimension(nsd) :: xi
      integer :: i,n1d,n
!
      select case(n)
!
!........................
!
      case(1)
!
!     constant
!
      dphi2dx=0.d0
!
!........................
!
      case(4)
!
!     bilinear
!
      n1d=2
!
      select case(i)
!
      case(1)
      dphi2dx=dphi(1,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(2)
      dphi2dx=dphi(2,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(3)
      dphi2dx=dphi(2,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(4)
      dphi2dx=dphi(1,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case default
!
      write(*,*) '4 node not defined'
      stop
!
      end select
!
!........................
!
      case(9)
!
!     biquadratic
!
      n1d=3
!
      select case(i)
!
      case(1)
      dphi2dx=dphi(1,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(2)
      dphi2dx=dphi(3,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(3)
      dphi2dx=dphi(3,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(4)
      dphi2dx=dphi(1,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(5)
      dphi2dx=dphi(2,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(6)
      dphi2dx=dphi(3,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(7)
      dphi2dx=dphi(2,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(8)
      dphi2dx=dphi(1,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(9)
      dphi2dx=dphi(2,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case default
!
      write(*,*) '5 node not defined', i
      stop
!
      end select
!
!........................
!
      case(16)
!
!     bicubic
!
      n1d=4
!
      select case(i)
!
      case(1)
      dphi2dx=dphi(1,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(2)
      dphi2dx=dphi(4,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(3)
      dphi2dx=dphi(4,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(4)
      dphi2dx=dphi(1,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(5)
      dphi2dx=dphi(2,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(6)
      dphi2dx=dphi(3,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(7)
      dphi2dx=dphi(4,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(8)
      dphi2dx=dphi(4,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(9)
      dphi2dx=dphi(3,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(10)
      dphi2dx=dphi(2,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(11)
      dphi2dx=dphi(1,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(12)
      dphi2dx=dphi(1,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(13)
      dphi2dx=dphi(2,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(14)
      dphi2dx=dphi(3,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(15)
      dphi2dx=dphi(3,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(16)
      dphi2dx=dphi(2,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case default
!
      write(*,*) '6 node not defined'
      stop
!
      end select
!
      case default
!
      write(*,*) 'degree not defined'
      stop
!
      end select
!
      end function
!
!-----------------------------------------------------------------------
!
      function dphi2dy(i,xi,n)
      use mgeometry, only : nsd
      implicit none
!
      real(8) :: dphi2dy
      real(8) :: phi,dphi
      real(8), dimension(nsd) :: xi
      integer :: i,n1d,n
!
      select case(n)
!
!........................
!
      case(1)
!
!     constant
!
      dphi2dy=0.d0
!
!........................
!
      case(4)
!
!     bilinear
!
      n1d=2
!
      select case(i)
!
      case(1)
      dphi2dy=phi(1,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(2)
      dphi2dy=phi(2,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(3)
      dphi2dy=phi(2,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(4)
      dphi2dy=phi(1,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case default
!
      write(*,*) '7 node not defined'
      stop
!
      end select
!
!........................
!
      case(9)
!
!     biquadratic
!
      n1d=3
!
      select case(i)
!
      case(1)
      dphi2dy=phi(1,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(2)
      dphi2dy=phi(3,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(3)
      dphi2dy=phi(3,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case(4)
      dphi2dy=phi(1,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case(5)
      dphi2dy=phi(2,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(6)
      dphi2dy=phi(3,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(7)
      dphi2dy=phi(2,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case(8)
      dphi2dy=phi(1,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(9)
      dphi2dy=phi(2,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case default
!
      write(*,*) '8 node not defined'
      stop
!
      end select
!
!........................
!
      case(16)
!
!     bicubic
!
      n1d=4
!
      select case(i)
!
      case(1)
      dphi2dy=phi(1,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(2)
      dphi2dy=phi(4,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(3)
      dphi2dy=phi(4,n1d,xi(1))*dphi(4,n1d,xi(2))
!
      case(4)
      dphi2dy=phi(1,n1d,xi(1))*dphi(4,n1d,xi(2))
!
      case(5)
      dphi2dy=phi(2,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(6)
      dphi2dy=phi(3,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(7)
      dphi2dy=phi(4,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(8)
      dphi2dy=phi(4,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case(9)
      dphi2dy=phi(3,n1d,xi(1))*dphi(4,n1d,xi(2))
!
      case(10)
      dphi2dy=phi(2,n1d,xi(1))*dphi(4,n1d,xi(2))
!
      case(11)
      dphi2dy=phi(1,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case(12)
      dphi2dy=phi(1,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(13)
      dphi2dy=phi(2,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(14)
      dphi2dy=phi(3,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(15)
      dphi2dy=phi(3,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case(16)
      dphi2dy=phi(2,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case default
!
      write(*,*) '9 node not defined'
      stop
!
      end select
!
      case default
!
      write(*,*) 'degree not defined'
      stop
!
      end select
!
      end function
!
!-----------------------------------------------------------------------
!
!     definida no elemento geometrico
!
      function psig(i,xx)
!      use mgeometry, only : nsd,nen, np
      use mgeometry
 
      implicit none
!
      real(8) :: psig
      real(8), dimension(nsd) :: xx
      integer :: i
!
      select case(i)
!
      case(1)
      psig=1.d0
!
      case(2)
      psig=xx(1)
!
      case(3)
      psig=xx(2)
!     
      case(4)
      psig = xx(1)*xx(2)
!     
      case(5)
      psig = xx(1)**2
!     
      case(6)
      psig = xx(2)**2
!
      case default
      write(*,*) 'basis function not defined'
      stop
      end select
!
      end function
!
!-----------------------------------------------------------------------
!
!     definida no elemento geometrico
!
      function dpsigx(i,xx)
!      use mgeometry, only : nsd,nen, np
      use mgeometry
 
      implicit none
!
      real(8) :: dpsigx
      real(8), dimension(nsd) :: xx
      integer :: i
!
      select case(i)
!
      case(1)
      dpsigx=0.d0
!
      case(2)
      dpsigx=1.d0
!
      case(3)
      dpsigx=0.d0
!     
      case(4)
      dpsigx = xx(2)
!     
      case(5)
      dpsigx = 2.d0*xx(1)
!     
      case(6)
      dpsigx = 0.d0
!
      case default
      write(*,*) 'basis function not defined'
      stop
      end select
!
      end function
!
!
!-----------------------------------------------------------------------
!
!     definida no elemento geometrico
!
      function dpsigy(i,xx)
!      use mgeometry, only : nsd,nen, np
      use mgeometry
 
      implicit none
!
      real(8) :: dpsigy
      real(8), dimension(nsd) :: xx
      integer :: i
!
      select case(i)
!
      case(1)
      dpsigy=0.d0
!
      case(2)
      dpsigy=0.d0
!
      case(3)
      dpsigy=1.d0
!     
      case(4)
      dpsigy=xx(2)
!     
      case(5)
      dpsigy=0.d0
!     
      case(6)
      dpsigy=2.d0*xx(2)
!
      case default
      write(*,*) 'basis function not defined'
      stop
      end select
!
      end function
!
!-----------------------------------------------------------------------
!
!     definida no elemento padrao 
!
      function thi(i, xi)
      use mgeometry, only : nsd,nen, np
      implicit none
!
      real(8) :: thi
      real(8), dimension(nsd) :: xi
      integer :: i,n1d
!
      select case(np)
!
      case(1)
!
      select case(i)
!
      case(1)
      thi = 1
!
      case default
      write(*,*) 'basis function not defined'
      stop
!
      end select
!
      case(3)
      select case(i)
      case(1)
      thi = 1
!
      case(2)
      thi = xi(1)
!
      case(3)
      thi = xi(2)
!
      case default
      write(*,*) 'basis function not defined'
      stop
!
      end select
!
      case default
!
      write(*,*) 'degree not defined'
      stop
!
      end select
      end function thi
!
!------------------------------------------------------
!
      subroutine jacobf(xe,xi,jf,detjf)
!
!     jacobian of the bilinear transformation
!
      use mgeometry
      implicit none
      real(8), dimension(nsd,nsd) :: jf
      real(8), dimension(nen,nsd) :: xe
      real(8), dimension(nsd)     :: xi
      real(8) :: detjf,eps=1.d-7!-6
!
      real(8) :: dphi2dx,dphi2dy
      integer :: i
!
      jf=0.d0
!
      do i=1,neng
!
      jf(1,1)=jf(1,1)+xe(i,1)*dphi2dx(i,xi,neng)
      jf(1,2)=jf(1,2)+xe(i,1)*dphi2dy(i,xi,neng)
!
      jf(2,1)=jf(2,1)+xe(i,2)*dphi2dx(i,xi,neng)
      jf(2,2)=jf(2,2)+xe(i,2)*dphi2dy(i,xi,neng)
!
      end do ! i
!

      detjf=jf(1,1)*jf(2,2)-jf(1,2)*jf(2,1)
!
      if(dabs(detjf).le.eps) then
      write(*,*) ' degenerated element'
      stop
      end if
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine velxy
!
      use mgeometry
      use marrays
      use mcoeficientes
!
      implicit none
!
      real(8) :: phi2d
      real(8), dimension(nsd)     :: vellocal, yy, xil
!
      integer :: i,j,k
!
      open(unit=70,file="vel.dat",status="unknown")
!
      do i=1,nelem
!
      yy(1) = (coord(el(i)%n(1),1) + coord(el(i)%n(2),1))/2.d0
      yy(2) = (coord(el(i)%n(1),2) + coord(el(i)%n(4),2))/2.d0
!
      vellocal = 0.d0
!
!     nodes of the element standart
! 
      do k=1,nsd
      xil(k)= 0.d0 ! ponto medio no elemento padrao
      end do ! k
!
!     velocity
!
      do k=1,nen
      vellocal(1) = vellocal(1) + phi2d(k,xil,nen)*bb(2*(el(i)%n(k))-1,1) !u
      vellocal(2) = vellocal(2) + phi2d(k,xil,nen)*bb(2*(el(i)%n(k)),1) !v
      end do ! k
!
      write(70   ,"(10(e25.15,2x))") yy(1), yy(2), vellocal(1), vellocal(2)
!
      end do !i 
!
      close(unit = 70)
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine prtgnuplot 
!
      use mgeometry
      use marrays
      use mtime, only: tt
!
      implicit none
!
      integer :: i,j,k
      real(8) :: zero=0.d0, pr, pe, h
      real(8), dimension(nsd) :: xx, yy
!
      open(unit = 10,file = 'elements.dat',status = 'unknown')
      open(unit = 20,file = 'solutionu.dat',status = 'unknown')
      open(unit = 30,file = 'solutionv.dat',status = 'unknown')
      open(unit = 40,file = 'solutionpres.dat',status = 'unknown')
      open(unit = 50,file = 'solutionprex.dat',status = 'unknown')
      open(unit = 60,file = 'solutionuex.dat',status = 'unknown')
      open(unit = 70,file = 'phi.dat',status = 'unknown')
!
      do i=1,nelem
!
      do j=1,4 !nen
!
      xx(1) = coord(el(i)%n(j),1)
      xx(2) = coord(el(i)%n(j),2)
!
      write(10,"(10(f25.15,2x))") coord(el(i)%n(j),1), coord(el(i)%n(j),2), zero
      write(20,"(10(f25.15,2x))") coord(el(i)%n(j),1), coord(el(i)%n(j),2), p(i,2*j-1) 
      write(30,"(10(f25.15,2x))") coord(el(i)%n(j),1), coord(el(i)%n(j),2), p(i,2*j) 
      write(40,"(10(f25.15,2x))") coord(el(i)%n(j),1), coord(el(i)%n(j),2), pressao(i,j) 
      write(50,"(10(f25.15,2x))") coord(el(i)%n(j),1), coord(el(i)%n(j),2), pr(xx) 
      write(60,"(10(f25.15,2x))") coord(el(i)%n(j),1), coord(el(i)%n(j),2), pe(xx,1)
      write(70,"(10(f25.15,2x))") coord(el(i)%n(j),1), coord(el(i)%n(j),2), h(xx,1.d0,tt)
!
      end do ! j
!
      xx(1) = coord(el(i)%n(1),1)
      xx(2) = coord(el(i)%n(1),2)
!
      write(10,"(10(f25.15,2x))") coord(el(i)%n(1),1), coord(el(i)%n(1),2), zero
      write(20,"(10(f25.15,2x))") coord(el(i)%n(1),1), coord(el(i)%n(1),2), p(i,1) 
      write(30,"(10(f25.15,2x))") coord(el(i)%n(1),1), coord(el(i)%n(1),2), p(i,2) 
      write(40,"(10(f25.15,2x))") coord(el(i)%n(1),1), coord(el(i)%n(1),2), pressao(i,1) 
      write(50,"(10(f25.15,2x))") coord(el(i)%n(1),1), coord(el(i)%n(1),2), pr(xx) 
      write(60,"(10(f25.15,2x))") coord(el(i)%n(1),1), coord(el(i)%n(1),2), pe(xx,1) 
      write(70,"(10(f25.15,2x))") coord(el(i)%n(1),1), coord(el(i)%n(1),2), h(xx,1.d0,tt) 
!
      write(10,*)
      write(10,*)
      write(20,*)
      write(20,*)
      write(30,*)
      write(30,*)
      write(40,*)
      write(40,*)
      write(50,*)
      write(50,*)
      write(60,*)
      write(60,*)
      write(70,*)
      write(70,*)
!
      end do ! i
!
      close(unit = 10)
      close(unit = 20)
      close(unit = 30)
      close(unit = 40)
      close(unit = 50)
      close(unit = 60)
      close(unit = 70)
!
      call velxy
!
      end subroutine
!
!!----------------------------------------------------------------------
!
      subroutine postprocessing
!
      use mgeometry
      use marrays
      use mcoeficientes
      use mtime, only:tt
!
      implicit none
!
      real(8), dimension(nsd)     :: xil,xx,xxg
      real(8), dimension(nen,nsd) :: xe,xer
!
      real(8) :: dphi2d,phi2d, thi,psig, h, hel
      integer :: i,j,k,n,ii,jj, nel, m
!
!      do i = 1,25
!      write(*,*) bb(i,1)
!      end do
!
!     coordinates of the nodes in the reference element
!
      select case(nen)
!
      case(4)
!
!     bilinear
!
      xer(1,1) = -1.d0
      xer(1,2) = -1.d0
      xer(2,1) =  1.d0
      xer(2,2) = -1.d0
      xer(3,1) =  1.d0
      xer(3,2) =  1.d0
      xer(4,1) = -1.d0
      xer(4,2) =  1.d0
!
      case(9)
!
!     biquadratic
!
      xer(1,1) = -1.d0
      xer(1,2) = -1.d0
      xer(2,1) =  1.d0
      xer(2,2) = -1.d0
      xer(3,1) =  1.d0
      xer(3,2) =  1.d0
      xer(4,1) = -1.d0
      xer(4,2) =  1.d0
!
      xer(5,1) =  0.d0
      xer(5,2) = -1.d0
      xer(6,1) =  1.d0
      xer(6,2) =  0.d0
      xer(7,1) =  0.d0
      xer(7,2) =  1.d0
      xer(8,1) = -1.d0
      xer(8,2) =  0.d0
      xer(9,1) =  0.d0
      xer(9,2) =  0.d0
!
      case default
!
      write(*,*) 'element degree not defined'
      stop
!
      end select
!
      p = 0.d0
      pressao = 0.d0
      cpres=0.d0
!
!     loop on elements
!
      do i=1,nelem
!
!     solution in the nodes of the element
!
      do j=1,nen
!
      do k=1,nsd
      xx (k)=xer(j,k)
      end do ! k
!
!     construction of solution in each element:
!
      do k = 1,nen
      p(i,2*j-1) = p(i,2*j-1) + phi2d(k,xx,nen)*bb(2*(el(i)%n(k))-1,1)
      p(i,2*j  ) = p(i,2*j  ) + phi2d(k,xx,nen)*bb(2*(el(i)%n(k))  ,1)
!      write(85,*) 2*(el(i)%n(k))-1, 2*(el(i)%n(k))
      end do ! k
      end do ! j
!
      end do ! i
!............................
!
!     encontra coeficientes que multiplicam funcoes de base da pressao
!
      do nel=1,nelem
!      cpres(:,nel) =-matmul(Api(nel,:,:)/lambda,p(nel,:))+Ag(nel,:,1)/lambda
      cpres(:,nel) =-matmul(Api(nel,:,:),p(nel,:))+Ag(nel,:,1)
!      write(171,*) nel,cpres(:,nel)
      end do !nel
!
      pressao=0.d0
!
!     encontrar pressao nos nos
!
      do i=1,nelem
!
!     solution in the nodes of the element
!
      do j=1,nen
      do k=1,nsd
      xe(j,k)=coord(el(i)%n(j),k)
      end do ! k
      end do ! j
!
      do j=1,nen
!
      do k=1,nsd
      xxg(k)= xe(j,k)
      end do ! k
!
!     construction of solution in each element:
!
      do k = 1,np
!
      pressao(i,j) = pressao(i,j) + psig(k,xxg)*cpres(k,i)
!      pressao(i,j) = pressao(i,j) + psig(k, xx)*cpres(k,i)
!
      end do ! k
!
!     construction of thickness matrix
!
      xx(1)=coord(el(i)%n(j),1)
      xx(2)=coord(el(i)%n(j),2)
      thick(i,j) = h(xx,pp,tt) ! h(coord(el(i)%n(j),1:nsd), pressao(i,j))
!
      end do ! j
!
      end do ! i
!
!      pressao=lambda*pressao

!     print '(e12.3)', pressao(4,:)
!     print '(e12.3)',  matmul(api(1,:,:),p(1,:))
!
      end subroutine
!

!----------------------------------------------------------------------
!
      subroutine vcheck
!
      use mgeometry
      use maverage
      use mcoeficientes
!
      implicit none
!
      real(8) :: v
!
      v=v0-(vol+volatot+volmtot)
!
      write(*,*) 'Volume inicial             :', v0
      write(*,*) 'Volume trocado pela aortica:', volatot
      write(*,*) 'Volume trocado pela mitral :', volmtot
      write(*,*) 'Volume vi - va - vm        :', v0-volatot-volmtot
      write(*,*) 'Volume atual               :', vol
      write(*,*) 'Check do balanco total     :', v,'=',dabs(v/vol)*1.d2,'%'
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine vvalvulas
!
      use mgeometry
      use mgauss
      use marrays
      use mcoeficientes
      use maverage
      use mtime, only:tt,dt
!
      implicit none
!
      integer :: nint=64,nint1d=8
      real(8), dimension(nen,nsd) :: xe
      real(8), dimension(4)       :: ee
      integer :: ii,jj,i,j,l,k,n1,n2,nel
      real(8) :: vola,volm,vole,xil1d,phi2d
      real(8) :: f4,f3,f7,flux
      real(8), dimension(nsd)     :: xil,xxl
!
      vola=0.d0
      volm=0.d0
!
!     aortica
!
      do jj=ny,ny
      do ii=1 ,nx
      nel=ii+nx*(jj-1)
!
      do j=1,nen
      do k=1,nsd
      xe(j,k) = coord(el(nel)%n(j),k)
      end do ! k
      end do ! j

      ee(1)=dsqrt((xe(1,1)-xe(4,1))**2.d0+(xe(1,2)-xe(4,2))**2.d0)
      ee(2)=dsqrt((xe(2,1)-xe(1,1))**2.d0+(xe(2,2)-xe(1,2))**2.d0)
      ee(3)=dsqrt((xe(3,1)-xe(2,1))**2.d0+(xe(3,2)-xe(2,2))**2.d0)
      ee(4)=dsqrt((xe(4,1)-xe(3,1))**2.d0+(xe(4,2)-xe(3,2))**2.d0)
!
      if(coord(el(nel)%n(7),1).ge.xai.and.coord(el(nel)%n(7),1).le.xaf) then
!
      vole=0.d0
!
      do l=1,nint1d
!
      xil1d = xig(l,nint1d)
      xil(1) = xil1d
      xil(2) = 1.d0
!
!     Real coordinates of the integration point
!
      xxl = 0.d0
!
      do n1=1,neng
      do n2=1,nsd
      xxl(n2) = xxl(n2) + phi2d(n1,xil,neng)*xe(n1,n2)
      end do ! n1
      end do ! n2
!
!     fluxo normal no ponto de integracao
!
      f4=p(nel, 8)
      f7=p(nel,14)
      f3=p(nel, 6)
      flux=f4*phi2d(4,xil,nen)+f7*phi2d(7,xil,nen)+f3*phi2d(3,xil,nen)
!
      vole=vole+flux*ee(4)/2.d0*wgl(l,nint1d)

      end do ! l
!
      vola=vola+vole
!
      end if
!
      end do ! ii
      end do ! jj
!
!     volume da aortica
!
      vola=vola*h0*dt
!
!     volume acumulado
!
      volatot=volatot+vola
!
!     mitral
!
      do jj=ny,ny
      do ii=1 ,nx
      nel=ii+nx*(jj-1)
!
      do j=1,nen
      do k=1,nsd
      xe(j,k) = coord(el(nel)%n(j),k)
      end do ! k
      end do ! j

      ee(1)=dsqrt((xe(1,1)-xe(4,1))**2.d0+(xe(1,2)-xe(4,2))**2.d0)
      ee(2)=dsqrt((xe(2,1)-xe(1,1))**2.d0+(xe(2,2)-xe(1,2))**2.d0)
      ee(3)=dsqrt((xe(3,1)-xe(2,1))**2.d0+(xe(3,2)-xe(2,2))**2.d0)
      ee(4)=dsqrt((xe(4,1)-xe(3,1))**2.d0+(xe(4,2)-xe(3,2))**2.d0)
!
      if(coord(el(nel)%n(7),1).ge.xmi.and.coord(el(nel)%n(7),1).le.xmf) then
!
      vole=0.d0
!
      do l=1,nint1d
!
      xil1d = xig(l,nint1d)
      xil(1) = xil1d
      xil(2) = 1.d0
!
!     Real coordinates of the integration point
!
      xxl = 0.d0
!
      do n1=1,neng
      do n2=1,nsd
      xxl(n2) = xxl(n2) + phi2d(n1,xil,neng)*xe(n1,n2)
      end do ! n1
      end do ! n2
!
!     fluxo normal no ponto de integracao
!
      f4=p(nel, 8)
      f7=p(nel,14)
      f3=p(nel, 6)
      flux=f4*phi2d(4,xil,nen)+f7*phi2d(7,xil,nen)+f3*phi2d(3,xil,nen)
!
      vole=vole+flux*ee(4)/2.d0*wgl(l,nint1d)

      end do ! l
!
      volm=volm+vole
!
      end if
!
      end do ! ii
      end do ! jj
!
!     volume da mitral
!
      volm=volm*h0*dt
!
!     volume acumulado
!
      volmtot=volmtot+volm
!

      write(64,'(f10.5,10(e15.8,1x))') tt,vola,volatot,volm,volmtot
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine volpres
!
      use mgeometry
      use mgauss
      use marrays
      use mcoeficientes
      use maverage
      use mtime, only:tt
!
      implicit none
!
      integer :: nint=64,nint1d=8
      real(8), dimension(nen,nsd) :: xe
      real(8), dimension(4)       :: ee
      real(8), dimension(nsd)       :: atoa
      real(8) :: h, area_el, area_tot, apres_ex, pr
      integer :: i,j,k, nel,l
      integer :: ii,jj,n1,n2,nk
      real(8), dimension(nsd,nsd) :: jf
      real(8) :: avol,vole,phi2d,detjf,apres,presme,ppm,presaor,presmit
      real(8), dimension(nsd)     :: xil,xxl
      real(8) :: f1,f2,f3,pp0,xil1d,pme
!
      avol =0.d0
      apres=0.d0
      apres_ex=0.d0
!
      vol  =0.d0
      presm=0.d0
!
      presaor=0.d0
      presmit=0.d0
!
      area_tot=(yint(2)-yint(1))*(xint(2)-xint(1))
!
      do nel = 1,nelem
!
      area_el=0.d0
!
      do j=1,nen
      do k=1,nsd
      xe(j,k) = coord(el(nel)%n(j),k)
      end do ! k
      end do ! j

      ee(1)=dsqrt((xe(1,1)-xe(4,1))**2.d0+(xe(1,2)-xe(4,2))**2.d0)
      ee(2)=dsqrt((xe(2,1)-xe(1,1))**2.d0+(xe(2,2)-xe(1,2))**2.d0)
!      ee(3)=dsqrt((xe(3,1)-xe(2,1))**2.d0+(xe(3,2)-xe(2,2))**2.d0)
!      ee(4)=dsqrt((xe(4,1)-xe(3,1))**2.d0+(xe(4,2)-xe(3,2))**2.d0)

      area_el = ee(1)*ee(2)
      avol    = avol  + area_el*h(coord(el(nel)%n(9),1:2),pp,tt)
      apres   = apres + area_el*pressao(nel,9)!*h(coord(el(nel)%n(9),1:2),pp,tt)
!      apres   = apres + area_el*pressao(nel,9) 
!
!     print*, tt, nel, p(nel,9)
!
!     integracao por quadratura gaussiana
!
      vole  =0.d0
      presme=0.d0
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k) = xi(l,nint,k)
      end do ! k
!
!     real coordinates of the integration point
!
      xxl=0.d0
!
      do j=1,neng
      do k=1,nsd
      xxl(k) = xxl(k) + phi2d(j,xil,neng)*xe(j,k)
      end do ! k
      end do ! j
!
!     pressao no ponto de integracao
!
      ppm=0.d0
      do j=1,neng
      ppm=ppm+phi2d(j,xil,neng)*pressao(nel,j)
      end do ! j
!
      call jacobf(xe,xil,jf,detjf)
!
      vole  =vole  +h(xxl,ppm,tt)*detjf*wg(l,nint)
!
      presme=presme+ppm*detjf*wg(l,nint)
!
      end do ! l
!
      vol  =vol  +vole
      presm=presm+presme
!
      end do !nel
!
      presm=presm/area_tot
!
!      apres = apres/avol
      apres = apres/area_tot!/0.036d0
!      if(tt<0.05)                            apres_ex = (186651.4d0*tt  + 1333.22d0)!
!      if(tt>=0.05.and.tt<0.3)                apres_ex = 10665.79d0+(tt-0.05d0)*42663.12d0+(tt-0.05d0)*(tt-0.175d0)*(-255978.56d0)
!      if(tt<=tiejec)                         apres_ex = (186651.4d0*tt  + 1333.22d0)!xx(2)
!      if(tt>tiejec.and.tt<=tirela)           apres_ex = 10665.79d0+(tt-0.05d0)*42663.12d0+(tt-0.05d0)*(tt-0.175d0)*(-255978.56d0)
                                                 !-55462.1d0*flowr/dt+ppa
                                                 ! 10665.79d0+(tt-0.05d0)*10665.8d0    
!      if(tt>tirela.and.tt<=tifill)           apres_ex = 13332.24d0-13332.24d0*(tt-0.3d0)/0.15d0 
!      if(tt>tifill.and.tt<=tcicle)           apres_ex = 1333.22d0*(tt-tifill)/(tcicle-tifill)
!      if(tt>tcicle) then
!      tt=tcicle
!      write(*,*) "tempo maior que um ciclo 4",tt
!      endif
!
      atoa(1)=0.d0
      atoa(2)=0.d0
      apres_ex=pr(atoa)
!
!     Pressao media na aortica
!
      nk=0
      do jj=ny,ny
      do ii=1 ,nx
      nel=ii+nx*(jj-1)
!
      do j=1,nen
      do k=1,nsd
      xe(j,k) = coord(el(nel)%n(j),k)
      end do ! k
      end do ! j

      ee(1)=dsqrt((xe(1,1)-xe(4,1))**2.d0+(xe(1,2)-xe(4,2))**2.d0)
      ee(2)=dsqrt((xe(2,1)-xe(1,1))**2.d0+(xe(2,2)-xe(1,2))**2.d0)
      ee(3)=dsqrt((xe(3,1)-xe(2,1))**2.d0+(xe(3,2)-xe(2,2))**2.d0)
      ee(4)=dsqrt((xe(4,1)-xe(3,1))**2.d0+(xe(4,2)-xe(3,2))**2.d0)
!
      if(coord(el(nel)%n(7),1).gt.(xai).and.coord(el(nel)%n(7),1).lt.(xaf)) then
!
      pme=0.d0
!
      do l=1,nint1d
!
      xil1d = xig(l,nint1d)
      xil(1) = xil1d
      xil(2) = 1.d0
!
!     Real coordinates of the integration point
!
      xxl = 0.d0
!
      do n1=1,neng
      do n2=1,nsd
      xxl(n2) = xxl(n2) + phi2d(n1,xil,neng)*xe(n1,n2)
      end do ! n1
      end do ! n2
!
!     fluxo normal no ponto de integracao
!
      f1 =pressao(nel,4)
      f2 =pressao(nel,7)
      f3 =pressao(nel,3)
      pp0=f1*phi2d(4,xil,nen)+f2*phi2d(7,xil,nen)+f3*phi2d(3,xil,nen)!
      pme=pme+pp0*ee(4)/2.d0*wgl(l,nint1d)

      end do ! l
!
      presaor=presaor+pme/ee(4)
      nk=nk+1
!
      end if
!
      end do ! ii
      end do ! jj
!
!     Pressao media
!
      presaor=presaor/nk
!
!     Pressao media na mitral
!
      nk=0.d0
      do jj=ny,ny
      do ii=1 ,nx
      nel=ii+nx*(jj-1)
!
      do j=1,nen
      do k=1,nsd
      xe(j,k) = coord(el(nel)%n(j),k)
      end do ! k
      end do ! j

      ee(1)=dsqrt((xe(1,1)-xe(4,1))**2.d0+(xe(1,2)-xe(4,2))**2.d0)
      ee(2)=dsqrt((xe(2,1)-xe(1,1))**2.d0+(xe(2,2)-xe(1,2))**2.d0)
      ee(3)=dsqrt((xe(3,1)-xe(2,1))**2.d0+(xe(3,2)-xe(2,2))**2.d0)
      ee(4)=dsqrt((xe(4,1)-xe(3,1))**2.d0+(xe(4,2)-xe(3,2))**2.d0)
!
      if(coord(el(nel)%n(7),1).ge.xmi.and.coord(el(nel)%n(7),1).le.xmf) then
!
      pme=0.d0
!
      do l=1,nint1d
!
      xil1d = xig(l,nint1d)
      xil(1) = xil1d
      xil(2) = 1.d0
!
!     Real coordinates of the integration point
!
      xxl = 0.d0
!
      do n1=1,neng
      do n2=1,nsd
      xxl(n2) = xxl(n2) + phi2d(n1,xil,neng)*xe(n1,n2)
      end do ! n1
      end do ! n2
!
!     fluxo normal no ponto de integracao
!
      f1 =pressao(nel,4)
      f2 =pressao(nel,7)
      f3 =pressao(nel,3)
      pp0=f1*phi2d(4,xil,nen)+f2*phi2d(7,xil,nen)+f3*phi2d(3,xil,nen)
!
      pme=pme+pp0*ee(4)/2.d0*wgl(l,nint1d)

      end do ! l
!
      presmit=presmit+pme/ee(4)
!
      nk=nk+1
!
      end if
!
      end do ! ii
      end do ! jj
!
!     Pressao media
!
      presmit=presmit/nk
!
!      write(63,'(f10.5,10(e15.8,1x))') tt, vol,avol, presm, apres, apres_ex
!
      write(63,'(f10.5,10(e15.8,1x))') tt, vol,apres_ex,presm,presaor,presmit
!      write(*,*) apres,presm
!
!stop
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine norms
!
      use mgeometry
      use mgauss
      use marrays
      use merro
!
      implicit none
!
      real(8), dimension(nsd)         :: xil,xx, yy
      real(8), dimension(nen,nsd)     :: xe
      real(8), dimension(nsd,nsd)     :: jf, invjf, tinvjf
      real(8) :: detjf
!
      real(8) :: pap             !pressao aproximada 
      real(8) :: pex             !pressao exata 
      real(8) :: elpres, el2pres !diferenca entre sol exata e aprox da pressao 
!
      real(8), dimension(2) :: pl, dplu, dplv    !sol exatas
      real(8), dimension(2) :: ph, dphu, dphv    !sol aproximadas
      real(8), dimension(2) :: elp, eldpu, eldpv !diferenca entre sol exata e aprox
      real(8) ::  el2p, el2dp      !norma dessas diferencas
!
      real(8) :: phi2d, dphi2dx, dphi2dy, thi,psig
      real(8) :: pe,dpex,dpey, pr
      integer :: i,j,k,l,ni,nj, nel
!
      integer :: nint
!
      nint=(8)**nsd
!
!     integrating
!
      el2p   = 0.d0
      el2dp  = 0.d0
      el2pres= 0.d0
!
!     loop on elements
!
      do i=1,nelem
!
!     coordinates of the element nodes
!
      do j=1,nen
      do k=1,nsd
      xe(j,k)=coord(el(i)%n(j),k)
      end do ! k
      end do ! j
!
!     loop on integration points
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k)=xi(l,nint,k)
      end do ! k
!
!     Jacobian of the transformation
!
      call jacobf(xe,xil,jf,detjf)
!
!     Inverse of the Jacobian
!
      invjf(1,1)= jf(2,2)
      invjf(1,2)=-jf(1,2)
      invjf(2,1)=-jf(2,1)
      invjf(2,2)= jf(1,1) 
!
      invjf = invjf/detjf
!
!     Transpose
!
      tinvjf = transpose(invjf)
!
!     real coordinates of the integration point
!
      xx=0.d0
!
      do ni=1,neng
      do nj=1,nsd
      xx(nj)=xx(nj)+phi2d(ni,xil,neng)*xe(ni,nj)
      end do ! nj
      end do ! ni
!
!     exact solutions at the integration points
!
      pex = pr(xx)
      pl(1) = pe(xx,1)
      pl(2) = pe(xx,2)
      dplu(1) = dpex(xx,1)
      dplu(2) = dpey(xx,1)
      dplv(1) = dpex(xx,2)
      dplv(2) = dpey(xx,2)
!
!     approximated solutions at the integration points
!
!      ph=ph+phip(k,xil)*el(i)%p(k,1)!
!
!     construct of solution in each element:
!
      ph   = 0.d0
      dphu = 0.d0
      dphv = 0.d0
      pap  = 0.d0
!
      do k = 1,nen
      ph(1) = ph(1) + phi2d(k,xil,nen)*bb(2*(el(i)%n(k))-1,1) !u
      ph(2) = ph(2) + phi2d(k,xil,nen)*bb(2*(el(i)%n(k)),1) !v
!
      dphu(1) = dphu(1) + ((tinvjf(1,1)*dphi2dx(k,xil,nen)+tinvjf(1,2)*dphi2dy(k,xil,nen)) &
                         )*bb(2*(el(i)%n(k))-1,1)
      dphu(2) = dphu(2) +  (tinvjf(2,1)*dphi2dx(k,xil,nen)+tinvjf(2,2)*dphi2dy(k,xil,nen) &
                         )*bb(2*(el(i)%n(k))-1,1)
!
      dphv(1) = dphv(1) + ((tinvjf(1,1)*dphi2dx(k,xil,nen)+tinvjf(1,2)*dphi2dy(k,xil,nen)) &
                         )*bb(2*(el(i)%n(k)),1)
      dphv(2) = dphv(2) + (tinvjf(2,1)*dphi2dx(k,xil,nen)+tinvjf(2,2)*dphi2dy(k,xil,nen) &
                         )*bb(2*(el(i)%n(k)),1)
!
      end do ! k

!      do nel=1,nelem
!      write(172,*) nel,cpres(:,nel)
!      end do !nel
      write(55,*) dplv(1), dphv(1), dplv(2), dphv(2)
      write(66,*) pl(1), ph(1), pl(2), ph(2)
!
      do k=1,np
      pap = pap + psig(k,xx)*cpres(k,i)
      end do !k
!
      elp    = dabs(ph-pl)
      eldpu  = dabs(dphu-dplu)     
      eldpv  = dabs(dphv-dplv)
      elpres = dabs(pex-pap)
!
!     l2 norm
!
      el2p  = el2p + (elp(1)**2+elp(2)**2)*wg(l,nint)*detjf
      el2pres  = el2pres + (elpres**2)*wg(l,nint)*detjf
!
!     grad norm
!
      el2dp = el2dp + (eldpu(1)**2 + eldpu(2)**2 + eldpv(1)**2 + eldpv(2)**2)*wg(l,nint)*detjf !tensor product : -->  integral de (grad(u)-grad(u_h)):(grad(u)-grad(u_h))
!
      end do ! l
!
      end do ! i
!
      el2p  = dsqrt(el2p)
      el2pres  = dsqrt(el2pres)
      el2dp = dsqrt(el2dp)
!
      open(unit = 100,file = 'norms.dat',status = 'unknown')
!
      write(*,"(15(e15.7))")  el2p, dlog(el2p), el2dp, dlog(el2dp), el2pres, dlog(el2pres)
      write(100,"(15(e15.7))") el2p, dlog(el2p), el2dp, dlog(el2dp), el2pres, dlog(el2pres)
!
      errol2p(ind,1)    = el2p
      errol2p(ind,2)    = dlog(el2p)
      errol2dp(ind,1)   = el2dp
      errol2dp(ind,2)   = dlog(el2dp)
      errol2pres(ind,1) = el2pres
      errol2pres(ind,2) = dlog(el2pres)
!
      end subroutine
!
!----------------------------------------------------------------------
!
!     sol exata : vetor u=(u,v)
!
      function pe(xx,i)
      use mgeometry, only : nsd, bc, xint, yint
      use mcoeficientes, only:mu
      use mtime, only: tt
      implicit none
!
      real(8), dimension(nsd) :: xx
      real(8) :: pe,pi, ter, um2
      integer :: i
      real(8) :: fluxtime
!
      pi = 4.d0*datan(1.d0)
!
      ter = xint(1) + (xint(2)-xint(1))/3.d0
      um2 = xint(1) + (xint(2)-xint(1))/2.d0
!
      select case(i)

!
      case(1)
!
      if (bc(1) == 1)      pe = (xx(2)-1.d0)*(xx(2)+1.d0)
!      if (bc(2) == 1)      pe = (xx(2)-yint(2))*(xx(2)-yint(1))/(2.d0*mu*(xint(1)-xint(2)))
!      if (bc(2) == 1)      pe = (1-xx(2)*xx(2))/(2.d0*mu) ! Inter [-1,1]
      if (bc(2) == 1)      pe = (1-xx(2))*xx(2)/(2.d0*mu) ! Inter [0,1]
      if (bc(3) == 1)      pe = -1.d0/(2.d0*pi*pi)*dsin(pi*xx(1))*dexp(xx(2)/2.d0)
      if (bc(4) == 1)      pe = pi*dcos(pi*xx(1))*dsin(pi*xx(2))
      if (bc(5) == 1) then
!      if (xx(2)==yint(1)) then
!      pe = 0.d0 !nao sei
!      else
      if(xx(1)==xint(1).or.xx(2)==yint(1).and.xx(1)<ter.or.xx(2)==yint(2)) then
      pe = 1.d0*fluxtime(tt)
      else
      pe = 0.d0 !nao usado
      end if
      end if
      if(bc(6)==1) then        
      pe = 0.d0
!      if(xx(1)==xint(1).and.xx(2)>0.5.and.xx(2).lt.yint(2)) pe=1.d0*fluxtime(tt)
      end if  
!
!...........
!
      case(2)
!
      if (bc(1) == 1)      pe = 0.d0
      if (bc(2) == 1)      pe = 0.d0
      if (bc(3) == 1)      pe = 1.d0/pi*dcos(pi*xx(1))*dexp(xx(2)/2.d0)
      if (bc(4) == 1)      pe = -pi*dsin(pi*xx(1))*dcos(pi*xx(2))
      if (bc(5) == 1) then
!      if (xx(2)==yint(1)) then
!      pe = 0.d0 !nao sei
!      else
      if(xx(1)==xint(1).or.xx(2)==yint(1).and.xx(1)<ter.or.xx(2)==yint(2)) then
      pe = 1.d0*fluxtime(tt)
      else
      pe = 0.d0 !nao usado
      end if
      end if
      if(bc(6)==1)         pe = 0.d0
!
!...........
!
      end select
!
      pe=0.d0
      end function
!
!----------------------------------------------------------------------
!
      function fluxtime(t)
!
      use mtime, only: tf
      implicit none
!
      real(8) :: fluxtime,t
      real(8) :: taux,eps,gamma
!  

      taux=tf!/3.d0!/2.d0!ti+0.2d0*tf
!      taux=0.1d0!
!
!      if(t.le.taux) then
!      fluxtime = 1.d0*(tt-ti)/(taux-ti)
!      fluxtime = 1.d0*dsin(pi/2.d0*(tt-ti)/(taux-ti))
!      else
!      fluxtime=1.d0
!      end if
!
!     fluxtime = 1-eps em tt=taux
!
      eps=1.d-3
      gamma=1.d0/taux*dlog((2.d0-eps)/eps)/2.d0
!
      fluxtime=dtanh(gamma*t)
      fluxtime=1.d0
!
      return
      end
!
!----------------------------------------------------------------------
!
!     exato: derivada de (u,v) em rel a x
!
      function dpex(xx,i)
      use mgeometry, only : nsd,bc
      implicit none
!
      real(8), dimension(nsd) :: xx
      real(8) :: dpex,pi
      integer :: i
!
      pi =4.d0*datan(1.d0)
!
      select case(i)
!
      case(1)
      if(bc(1)==1)      dpex = 0.d0 
      if(bc(2)==1)      dpex = 0.d0 
      if(bc(3)==1)      dpex = -dcos(pi*xx(1))*dexp(xx(2)/2)/(2*pi) !Correa
      if(bc(4)==1)      dpex = -pi*pi*dsin(pi*xx(1))*dsin(pi*xx(2))
      if(bc(5)==1)      dpex = 0.d0 !nao sei
      if(bc(6)==1)      dpex = 0.d0 !nao sei
!
      case(2)
!
      if(bc(1)==1)      dpex = 0.d0 
      if(bc(2)==1)      dpex = 0.d0 
      if(bc(3)==1)      dpex = -dsin(pi*xx(1))*dexp(xx(2)/2) !Correa
      if(bc(4)==1)      dpex = -pi*pi*dcos(pi*xx(1))*dcos(pi*xx(2))
      if(bc(5)==1)      dpex = 0.d0 !nao sei
      if(bc(6)==1)      dpex = 0.d0 !nao sei
!
      end select
!
      end function
!
!----------------------------------------------------------------------
!
!     exato: derivada de (u,v) em rel a y
!
      function dpey(xx,i)
      use mgeometry, only : nsd,bc, xint, yint
      use mcoeficientes, only:mu
      implicit none
!
      real(8), dimension(nsd) :: xx
      real(8) :: dpey,pi
      integer :: i
!
      pi =4.d0*datan(1.d0)
	
      select case(i)
!
!...........
!
      case(1)
!
      if(bc(1)==1)      dpey = 2.d0*xx(2) 
!      if(bc(2)==1)      dpey = 1.d0/(2.d0*mu*(xint(2)-xint(1)))*(-yint(2)-yint(1)+2*xx(2))
      if(bc(2)==1)      dpey = (1.d0-2.d0*xx(2))/(2.d0*mu)
      if(bc(3)==1)      dpey = -dsin(pi*xx(1))*dexp(xx(2)/2)/(4*pi*pi)
      if(bc(4)==1)      dpey = pi*pi*dcos(pi*xx(1))*dcos(pi*xx(2))
      if(bc(5)==1)      dpey = 0.d0 !nao sei
      if(bc(6)==1)      dpey = 0.d0 !nao sei
!
!............
!
      case(2)
!
      if(bc(1)==1)      dpey = 0.d0 
      if(bc(2)==1)      dpey = 0.d0 
      if(bc(3)==1)      dpey = dcos(pi*xx(1))*dexp(xx(2)/2)/(2*pi)  !escrever
      if(bc(4)==1)      dpey = pi*pi*dsin(pi*xx(1))*dsin(pi*xx(2))
      if(bc(5)==1)      dpey = 0.d0 !nao sei
      if(bc(6)==1)      dpey = 0.d0 !nao sei
!
      end select
!
      end function
!
!----------------------------------------------------------------------
!
!     exato: distribuicao de pressao
!
      function pr(xx)
      use mgeometry
      use mcoeficientes
      use mtime, only: tt, tf, dt
      use marrays, only:ppa, flowr
!
      implicit none
      real(8), dimension(nsd) :: xx
      real(8), dimension(nsd) :: pr0, pr1, pr2, pr3, pr4, pr5, pr6
      real(8) :: pr,pi
      real(8) :: dpex, dpey, fluxtime
      real(8) :: interpol1, interpol2
!
      pi = 4.d0*datan(1.d0)
!
      if(bc(1)==1)                           pr = dsin(pi*xx(1))*(xx(2)-1.d0)*(xx(2)+1.d0)
!      if(bc(2)==1)                           pr = (xx(1)-xint(2))/(xint(1)-xint(2))
      if(bc(2)==1)                           pr = 1.d0-xx(1) ! Inter = [0,1]
      if(bc(3)==1)                           pr = -1.d0/pi*dcos(pi*xx(1))*dexp(xx(2)/2.d0) ! Correa 
      if(bc(4)==1)                           pr = (1.d0+2.d0*pi*pi)*dsin(pi*xx(1))*dsin(pi*xx(2))
      if(bc(5)==1)                           pr = 0.d0 !tracao sendo nula
!
!
      pr0(1)=t0
      pr0(2)=p0
      pr1(1)=t1
      pr1(2)=p1
      pr2(1)=t2
      pr2(2)=p2
      pr3(1)=t3
      pr3(2)=p3
      pr4(1)=t4
      pr4(2)=p4
      pr5(1)=t5
      pr5(2)=p5
      pr6(1)=t6
      pr6(2)=p6
!
!
      if(bc(6)==1)   then
!
      if(tt<=t1)             pr = interpol1(tt,pr0,pr1)
      if(tt>t1.and.tt<=t3)   pr = interpol2(tt,pr1,pr2,pr3)     
      if(tt>t3.and.tt<=t4)   pr = interpol1(tt,pr3,pr4)
      if(tt>t4.and.tt<=t6)   pr = interpol1(tt,pr4,pr6)!interpol2(tt,pr4,pr5,pr6)
      if(tt>t6) then
      tt=t6
      write(*,*) "tempo maior que um ciclo 5",tt
      endif
      endif
!
      end function
!
!----------------------------------------------------------------------
!
!     Interpolation with lagrangian polynomial order 1
!
      function interpol1(t,v0,v1)
!
      implicit none
      real(8), dimension(2) :: v0,v1
      real(8) :: t, interpol1
!
      interpol1= (v1(2)*(t-v0(1))+v0(2)*(v1(1)-t))/(v1(1)-v0(1))
!
      end function
!
!----------------------------------------------------------------------
!
!     Derivative in t of the interpolation with lagrangian polynomial order 1
!
      function interpol1dt(v0,v1)
!
      implicit none
      real(8), dimension(2) :: v0,v1
      real(8) :: t, interpol1dt
!
      interpol1dt= (v1(2)-v0(2))/(v1(1)-v0(1))
!
      end function
!----------------------------------------------------------------------
!
!     Interpolation with lagrangian polynomial order 1
!
      function interpol2(t,v0,v1,v2)
!
      implicit none
      real(8), dimension(2) :: v0,v1,v2
      real(8) :: t, interpol2
!
      interpol2= v0(2)*(t-v1(1))*(t-v2(1))/(v0(1)-v1(1))/(v0(1)-v2(1)) + &
               & v1(2)*(t-v0(1))*(t-v2(1))/(v1(1)-v0(1))/(v1(1)-v2(1)) + &
               & v2(2)*(t-v0(1))*(t-v1(1))/(v2(1)-v0(1))/(v2(1)-v1(1))
!
      end function
!
!---------------------------------------------------------------------------
!
!     subrotina de inversao                                              
!
      subroutine invmb(am,m)                                             
!
      implicit real*8(a-h,o-z)                                           
      dimension ipi(m),ind(m,2),piv(m),dis(m,1),am(m,m)    
!                 
      ncoln=0                                                            
      det=1.0                                                            
      do 20 j=1,m                                                        
   20 ipi(j)=0                                                           
      do 550 i=1,m                                                       
      amax=0.0                                                           
      do 105 j=1,m                                                       
      if(ipi(j)-1)60,105,60                                              
   60 do 100 k=1,m                                                       
      if(ipi(k)-1) 80,100,740                                            
   80 if(dabs(amax)-dabs(am(j,k)))85,100,100                               
   85 irow=j                                                             
      ico=k                                                              
      amax=am(j,k)                                                       
  100 continue                                                           
  105 continue                                                           
      ipi(ico)=ipi(ico)+1                                                
      if(irow-ico)140,260,140                                            
  140 det=-det                                                           
      do 200 l=1,m                                                       
      swap=am(irow,l)                                                    
      am(irow,l)=am(ico,l)                                               
  200 am(ico,l)=swap                                                     
      if(ncoln) 260,260,210                                              
  210 do 250 l=1,ncoln                                                   
      swap=dis(irow,l)                                                   
      dis(irow,l)=dis(ico,l)                                             
  250 dis(ico,l)=swap                                                    
  260 ind(i,1)=irow                                                      
      piv(i)=am(ico,ico)                                                 
      ind(i,2)=ico                                                       
      det=det*piv(i)                                                     
      am(ico,ico)=1.0                                                    
      do 350 l=1,m                                                       
  350 am(ico,l)=am(ico,l)/piv(i)                                         
      if(ncoln) 380,380,360                                              
  360 do 370 l=1,ncoln                                                   
  370 dis(ico,l)=dis(ico,l)/piv(i)                                       
  380 do 550 lz=1,m                                                      
      if(lz-ico)400,550,400                                              
  400 t=am(lz,ico)                                                       
      am(lz,ico)=0.0                                                     
      do 450 l=1,m                                                       
  450 am(lz,l)=am(lz,l)-am(ico,l)*t                                      
      if(ncoln)550,550,460                                               
  460 do 500 l=1,ncoln                                                   
  500 dis(lz,l)=dis(lz,l)-dis(ico,l)*t                                   
  550 continue                                                           
      do 710 i=1,m                                                       
      l=m+1-i                                                            
      if(ind(l,1)-ind(l,2))630,710,630                                   
  630 jrow=ind(l,1)                                                      
      jco=ind(l,2)                                                       
      do 705 k=1,m                                                       
      swap=am(k,jrow)                                                    
      am(k,jrow)=am(k,jco)                                               
      am(k,jco)=swap                                                     
  705 continue                                                           
  710 continue                                                           
  740 continue                                                           
      return                                                             
      end
!--------------------------------------------------------------
