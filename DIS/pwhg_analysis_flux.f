!  The next subroutines, open some histograms and prepare them 
!      to receive data 
!  You can substitute these  with your favourite ones
!  init   :  opens the histograms
!  topout :  closes them
!  pwhgfill  :  fills the histograms with data

      subroutine init_hist
         implicit none
         include  'LesHouches.h'
         include 'pwhg_math.h'
         include 'nlegborn.h'
         include 'pwhg_kn.h'
         call inihists
         
         call bookupeqbins('sigtot',1d0,0d0,1d0)
         call bookupeqbins('Q',  0.5d0, 2d0, 20d0)

         call bookupeqbins('Q2', (1000d0-25d0)/20d0, 25d0, 1000d0)
         call bookupeqbins('y',  (0.95d0-0.04d0)/20d0 , 0.04d0, 0.95d0)
         call bookupeqbins('x',  0.05d0 , 0d0, 1d0)
         call bookupeqbins('xl',  0.05d0 , 0d0, 1d0)   
         call bookupeqbins('ptl', 2d0, 0d0, 100d0)
         call bookupeqbins('etal',0.1d0,0d0,5d0)
         call bookupeqbins('ptj_overQ_breit', 0.01d0, 0d0, 0.5d0)

      end
       
      subroutine analysis(dsig0)
         implicit none
         real * 8 dsig0
         include 'nlegborn.h'
         include 'hepevt.h'
         include 'pwhg_math.h' 
         include 'LesHouches.h'
         include 'pwhg_kn.h'
         include 'pwhg_weights.h'
         ! include 'pwhg_rad.h'
         ! include 'pwhg_rwl.h'
         logical ini
         data ini/.true./
         save ini
         real * 8 dsig(weights_max)
   
         integer   maxtrack,maxjet
         integer i,njets
         parameter (maxtrack=2048,maxjet=128)
         ! we need to tell to this analysis file which program is running it
         character * 6 WHCPRG
         common/cWHCPRG/WHCPRG
         ! data WHCPRG/'NLO   '/
         integer nlep, npartons,nunident                         !lepton, quark, lepton initialstate, parton initialstate 
         real * 8 plep(4,maxtrack), plephard(4),pquark(4,maxtrack),plis(4),ppis(4),q(4)              !momenta of lepton, quark, lepton initialstate, parton is, momentum transfere
         real * 8 El_hardest
         real * 8 plab(0:3,maxtrack), pjets(0:3,maxjet)
         real * 8 y,x,Q2, xl               !xs. inv mass of incoming parton-electron, momentum transfer variable Q2=-q^2
         real * 8 ptj(maxjet),yj(maxjet),phij(maxjet), eta1, ptjetmin, absEtaMax
         real * 8 ptl, etal
         real * 8, external :: phepdot, eta, kt2
         real * 8 refin(1:4), refout(1:4)
         real * 8 alpha, beta, ptbreit
         real * 8, save :: Eproton, Elepton, sbeams
         logical, save :: fixed_target
         real *8, external :: powheginput
         
       
        
         dsig=0d0
         call multi_plot_setup(dsig0,dsig, weights_max)
   
         if(ini) then         
            if(whcprg.eq.'NLO') then
               write(*,*) '       NLO analysis'
            elseif(WHCPRG.eq.'LHE   ') then
               write(*,*) '       LHE analysis'
            elseif(WHCPRG.eq.'HERWIG') then
               write (*,*) '           HERWIG ANALYSIS            '
            elseif(WHCPRG.eq.'PYTHIA') then
               write (*,*) '           PYTHIA ANALYSIS            '
            endif
            ini=.false.

            Elepton = powheginput("ebeam1")
            Eproton = powheginput("ebeam2")
            fixed_target = (powheginput("#fixed_target") .eq. 1d0)

            if(fixed_target) then
               sbeams = 2*Elepton*Eproton + Eproton**2
            else
               sbeams = 4 * Elepton * Eproton
            endif


            if(whcprg.eq.'NLO') then
               if( abs( 4d0 * ebmup(1) * ebmup(2) - sbeams) .gt. 1d-4)
     $              then
                  write(*,*) "inconsistency in the calculation of ",
     $                 "the partonic center-of-mass energy in the ",
     $                 "analysis, aborting ..."
                  stop
               endif
            endif
            
         endif
   
         nlep = 0
         npartons = 0
         pquark = 0d0
         plep = 0d0
         ppis = 0d0
         plis = 0d0
         pquark = 0d0
         nunident = 0
         if(whcprg.eq.'NLO   '.or.whcprg.eq.'LHE   '
     $        .or.whcprg.eq.'PYTHIA') then
            do i=1,nhep
               !if(idhep(i).eq.11) print*, phep(1:4,i), eta(phep(:,i)), sqrt(kt2(phep(:,i)))
   !     Initial states
               if(isthep(i).eq.-1.or.isthep(i).eq.21) then
                  if(abs(idhep(i)).le.16 .and.abs(idhep(i)).ge.11 ) then
                     plis(1:4) = phep(1:4,i)
                  else if(abs(idhep(i)).le.5.or.idhep(i).eq.21) then
                     ppis(1:4) = phep(1:4,i)
                  else
                     stop 'Wrong initial state'
                  endif
   !     Final states
               else if(isthep(i).eq.1) then
                  if(abs(idhep(i)).ge.11 .and. abs(idhep(i)).le.16) then
                     nlep = nlep + 1
                     plep(1:4,nlep) = phep(1:4,i)
                  elseif (abs(idhep(i)) <= 9 .or. abs(idhep(i)) == 21 .or. abs(idhep(i)) > 100 ) then 
                     npartons = npartons + 1
                     pquark(1:4,npartons) = phep(1:4,i)
                  else
                     nunident = nunident + 1
                     ! print*, 'idhep(i)', idhep(i)
                  endif
               endif
            enddo
         else
            print*, 'Wrong program', whcprg
            stop
         endif

         if (nunident>0) print*, 'nunident', nunident
   
         if(nlep<1) stop 'Wrong number of leptons'
   
         El_hardest = 0d0
         do i = 1, nlep
            if (plep(4,i) > El_hardest ) then
               El_hardest = plep(4,i)
               plephard(1:4) = plep(1:4,i) 
            endif
         end do
   
   
         q(:) = plis(:) - plephard(:)

         xl = plis(4)/ebmup(1)
         
         Q2 = phepdot(q,q)
         Q2 = -Q2

!     Q^2 = xl xdis y S
!     x = Q^2/(S*y*xdis)
!     y = (q.P)/(lin.P)
         
         y = phepdot(ppis,q) / phepdot(ppis,plis)      
         x = Q2 / (sbeams * y * xl)
    
         
         call filld('sigtot',0.5d0,dsig)   
         call filld('Q2', Q2, dsig)
         call filld('Q', sqrt(Q2), dsig)
         call filld('x', x, dsig)
         call filld('y', y, dsig)
         call filld('xl', xl, dsig)

         ptl  = sqrt(kt2(plephard(:)))
         call getrapidity(plephard(:),etal)
         call filld('ptl', ptl, dsig)
         call filld('etal', etal, dsig)


         if(npartons == 2) then
            refin  = x *ebmup(2) * (/ 0d0, 0d0, -1d0, 1d0/)
            if(fixed_target) refin = refin/2d0 ! The fake massless beam has E = mp/2d0, hence the factor 2
            if(plis(3) < 0) refin(3)=-refin(3) ! Flip the sign of z if necessary
            refout = q(:) + refin
            
            alpha = 2d0 * phepdot(pquark(:, 1), refout)/Q2
            beta  = 2d0 * phepdot(pquark(:, 1), refin)/Q2
            
            ptbreit = sqrt(alpha*beta)
         else
            ptbreit = 0d0
         endif
         call filld('ptj_overQ_breit', ptbreit, dsig)

         
         

         
c$$$         print *, "Check event!"
c$$$         print *, "npartons = ", npartons
c$$$         print *, "Q2 = ", Q2, 2* phepdot(refin, refout)
c$$$         print *, "refout.m2 = ", phepdot(refout,refout) 
  
      end
  
  
  
      function phepDot(p_A,p_B)
         implicit none
         real * 8  phepDot
         real * 8  p_A(4),p_B(4)
         phepDot=p_A(4)*p_B(4)-p_A(1)*p_B(1)-p_A(2)*p_B(2)-p_A(3)*p_B(3)
      end
  
      function kt2(p)
         implicit none
         real * 8 kt2, p(1:4)
   
         kt2 = p(1)**2 + p(2)**2
      end
  
      function eta(p)
         implicit none
         real * 8 eta, p(0:3), normp, norm
   
         normp = norm(p)
         if(normp.gt.p(3)) then
            eta = 0.5d0 * log((normp + p(3)) / (normp - p(3)))
         else
            eta = sign(1d100,p(3)) 
         endif
      end
  
      function norm(p)
         implicit none
         real * 8 norm, p(0:3)
   
         norm = p(1)**2 + p(2)**2 + p(3)**2
         norm = sqrt(max(0d0,norm))
         end
   
         function getrapidity0(p)
         implicit none
         real * 8 p(0:3),getrapidity0
         getrapidity0=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end
  
      subroutine getrapidity(p,y)
         implicit none
         real * 8 p(4),y
         y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      end
  
      function azi(p)
         implicit none
         real * 8 pi
         parameter(pi = 3.141592653589793D0)
         real * 8 azi,p(0:3)
         azi = atan(p(2)/p(1))
         if (p(1).lt.0d0) then
            if (azi.gt.0d0) then               
               azi = azi - pi
            else
               azi = azi + pi
            endif
         endif    
      end
        
      ! Builds pp like jets. Should take as argument all final state
      ! partons/hadrons. Returns the jets pt ordered. 
      subroutine buildjets(n,pin,pj,njets,ptj,yj,phij, ptjetmin, absEtaMax)
         implicit none
         integer n
         double precision pin(0:3,n) 
         integer maxtrack,maxjet
         parameter (maxtrack=2048,maxjet=128)
   
   !     Output
         double precision pj(0:3,maxjet)
   !     Internal
         integer mu, njets, njets_all, ntracks, ijet, j
         integer jetvec(maxtrack)
         double precision pjet(4,maxjet), pjet_all(4,maxjet)
         double precision ptrack(4,maxtrack)
         double precision ptj(maxjet),yj(maxjet),phij(maxjet), ptjetmin, absEtaMax
         double precision ptj_all(maxjet),yj_all(maxjet),phi_all(maxjet)
         double precision R, ptmin_fastkt, palg
         double precision, external :: azi, eta
   
   
         ptrack = 0d0
         jetvec = 0
         pjet = 0d0
         pjet_all = 0d0
         njets=0
         njets_all = 0
         ntracks = n
         
         ptrack(4,1:n)=pin(0,1:n)
         ptrack(1:3,1:n)=pin(1:3,1:n)
         
         R = 0.8d0
         ptmin_fastkt = 0d0
         palg = -1d0
   !      palg = 0d0
   !      palg = -1d0
c     -1 is anti_kt, +1 is kt
  
         call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin_fastkt,pjet,njets,jetvec)

         if(njets.gt.maxjet) then
            print*, 'njets out of bounds!!', njets, maxjet
            stop
         endif
         
         j=0
         pjet_all=pjet
         njets_all=njets
         pjet=0d0
         njets=0d0
         ptj = 0d0
         yj = 0d0
         
         do ijet=1,njets_all
            ptj_all(ijet) = sqrt(pjet_all(1,ijet)**2 + pjet_all(2,ijet)**2)
            call getrapidity(pjet_all(:,ijet),yj_all(ijet))
            phi_all(ijet) = azi(pjet_all(:,ijet))
  
          if(ptj_all(ijet).gt.ptjetmin.and.
     $        abs(eta(cshift(pjet_all(:,ijet),-1))).lt.absEtaMax) then
      !     if(ptj_all(ijet).gt.ptalljetmin.and.
      ! $        abs(yj_all(ijet)).lt.yjetmax) then
               j=j+1
               pjet(:,j)=pjet_all(:,ijet)
               ptj(j) = ptj_all(ijet)
               yj(j) = yj_all(ijet)
               phij(j) = phi_all(ijet)
          endif
         enddo
         njets=j
   
         pj = 0d0
   
         do ijet=1,njets
            do mu=1,3
               pj(mu,ijet)=pjet(mu,ijet)
            enddo
            pj(0,ijet)=pjet(4,ijet)
         enddo
              
      end
  
  
