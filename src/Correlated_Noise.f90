
module global_vars 
	implicit none
	! Global run configuration, read once from input/input.dat.
	! Shared by main and helper routines that `use global_vars`.
	integer :: iseed, isim !seed
	integer :: idir, idir_read, nsample, ntau !idir=1 for directed
	real(8) :: DD, dt !noise, dt
	real(8) :: t0, tend, prob !start, end, gij, prob
	real(8) :: W_mean, W_sig !mean and sig of weighted random graph
	real(8) :: r0

	real(8), parameter :: pi=4d0*atan(1d0)

	! === 動態網路設定 ===
	integer :: NN                      ! 動態節點數（取代 parameter N）
	integer :: use_file_network        ! 0 = 程式生成網路, 1 = 從檔案讀取 Aij
	character(len=256) :: network_file ! 網路檔案路徑
	
	contains

	subroutine read_variables_from_file()
	
	  open(unit=10, file='input/input.dat', status='old')
	  read(10, *) dt, iseed, isim ! simulation settings
	  read(10, *) t0, tend, nsample, ntau  ! simulation settings
      ! t0, tend, icon, ntau 分別是模擬的開始時間、結束時間、每隔多少步取樣一次、<x(t+tau)x(t)>的tau值
	  read(10, *) idir, prob, r0, W_mean, W_sig, DD ! network settings
      ! idir=1 for directed, prob是連結機率，r0是耦合強度的縮放參數，W_mean和W_sig是加權隨機圖的平均和標準差，dd是噪聲強度
	  read(10, *) use_file_network, idir_read, NN ! 新增：0=生成, 1=讀檔；以及節點數
	  if (use_file_network .eq. 1) then
	    read(10, '(A)') network_file   ! 新增：網路檔案路徑
	    network_file = adjustl(network_file) ! adjustl的作用是把字串開頭的空白去掉，移到字串尾端。
	  endif
	  close(10)
	end subroutine read_variables_from_file
	
	end module global_vars

    program Correlated_Noise
    use global_vars
    implicit none
    integer, parameter :: MAXT_Ktau = 5
    integer :: i, j, k, ii, jj, ia, ib
    integer :: i_avg, ITIME, Icount, jcount, IT, JTEM
    integer :: N, NSTEP, nedge1
    integer :: info
    real(8) :: tt, dum, wij, mean_degree
    real(8) :: time1, time2

    integer, allocatable :: G0(:,:), indegree(:), outdegree(:)
    real(8), allocatable :: x1(:), f1(:), gau(:), noise_increment(:)
    real(8), allocatable :: r(:), noise(:,:), noise_corr(:,:), noise_diag(:,:), noise_cholesky(:,:)
    real(8), allocatable :: G(:,:), Q(:,:), Lyapunov_check(:,:)
	real(8), allocatable :: inweightdegree(:), outweightdegree(:)

    real(8), allocatable :: avex(:), avef(:), Ko(:,:), Ktau(:,:,:), KtauBUF(:,:)

    real(8), allocatable :: QschurSol(:,:), QschurWORK(:)
    real(8), allocatable :: QschurT(:,:), QschurVS(:,:), QschurWR(:), QschurWI(:)
	logical, allocatable :: QschurBWORK(:)
    complex(8), allocatable :: KoT(:,:) 
    character(len=256) :: header_str  ! 讀網路檔標頭用
    real(8), external :: RAN3, GAUDEV

    call read_variables_from_file()
	N = NN

    allocate(G0(N,N), indegree(N), outdegree(N))
    allocate(x1(N), f1(N), gau(N))
    allocate(r(N), noise(N,N), noise_corr(N,N), noise_diag(N,N), noise_cholesky(N,N), noise_increment(N))
    allocate(G(N,N), Q(N,N), inweightdegree(N), outweightdegree(N))
    allocate(avex(N), avef(N), Ko(N,N), Ktau(N,N,MAXT_Ktau), KtauBUF(N,0:MAXT_Ktau-1))
    allocate(QschurT(N,N), QschurVS(N,N), QschurSol(N,N))
    allocate(QschurWR(N), QschurWI(N), QschurWORK(max(1,8*N)), QschurBWORK(N))
    allocate(KoT(N,N))

    OPEN(UNIT=11, FILE='config/parameter.csv', STATUS='REPLACE')
    OPEN(UNIT=12, FILE='config/noise.csv', STATUS='REPLACE')
    OPEN(UNIT=13, FILE='config/G_connect.csv', STATUS='REPLACE')

    OPEN(UNIT=14, FILE='output/average/x.dat', STATUS='REPLACE')
    OPEN(UNIT=15, FILE='output/average/F.dat', STATUS='REPLACE')
    OPEN(UNIT=16, FILE='output/average/Ko.dat', STATUS='REPLACE')
    OPEN(UNIT=17, FILE='output/average/Ktau.dat', STATUS='REPLACE')

    OPEN(UNIT=18, FILE='output/solution/Lyapunov_check.dat', STATUS='unknown')
    OPEN(UNIT=19, FILE='output/solution/Koii.dat', STATUS='unknown')
	OPEN(UNIT=20, FILE='output/solution/Koij.dat', STATUS='unknown')

    mean_degree = prob*N

    write(6,*) repeat("=", 40)
	write(6,*) 'number of node = ', N
    write(6,*) 'mean degree = ', mean_degree
	write(6,*) repeat("=", 40)

    iseed = -iseed; dum = RAN3(ISEED)
	nstep = INT((tend-t0)/dt)	

    G = 0d0; G0 = 0

    ! Convention: G(i,j) affects node i from node j.
	if (use_file_network .eq. 1) then
	! === 從 Aij 檔讀取網路 ===
		print*, 'Reading network from file: ', trim(network_file)
		! trim() 是去掉字串尾端多餘的空白
		open(unit=999, file=trim(network_file), status='old', action='read')
		read(999, *, err=200, end=200) header_str  ! 嘗試跳過標頭行
		do
			read(999, *, end=100) jj, ii, wij ! 注意這裡讀取順序是 jj, ii，對應到 G(ii, jj)
			if (ii .ge. 0 .and. ii .le. N-1 .and. jj .ge. 0 .and. jj .le. N-1) then
				G(ii+1, jj+1) = wij
				if (wij .gt. 0d0) G0(ii+1, jj+1) = 1
			endif
		end do
100		continue
		close(999)
		goto 201
200		continue
		! 若標頭行讀取失敗，重新開啟從頭讀（無標頭格式）
		close(999)
		open(unit=999, file=trim(network_file), status='old', action='read')
		do
			read(999, *, end=201) ii, jj, wij
			if (ii .ge. 1 .and. ii .le. N .and. jj .ge. 1 .and. jj .le. N) then
				G(ii, jj) = wij
				if (wij .gt. 0d0) G0(ii, jj) = 1
			endif
		end do
201		continue
		close(999)
		if (idir_read .eq. 0) then
		! 若讀取的網路是無向的，則對稱填充 G 和 G0
		G = G + transpose(G)
		G0 = G0 + transpose(G0)
		endif
	else
	! === 程式生成網路 ===
		if (idir .eq. 0) print*, 'Undirected ER random graph'
		if (idir .eq. 1) print*, 'Directed ER random graph'	 		
		do i = 1, N-1; do j = i+1, N
		   	G0(i,j)=INT(1.d0+prob-ran3(iseed))
		   	if(idir.eq.0) G0(j,i)=G0(i,j)
		   	if(idir.eq.1) G0(j,i)=INT(1.d0+prob-ran3(iseed))	      
		enddo; enddo	

		if (idir .eq. 0) then
			do i = 1, n-1; do j = i+1, n
			! G(i,j)=gau_mean*G0(i,j)
			G(i,j)=GAUDEV(W_mean, W_sig, iseed)*G0(i,j)
			G(j,i)=G(i,j)
			enddo; enddo
		else
			do i = 1, n-1; do j = i+1, n
			! G(i,j)=gau_mean*G0(i,j) 
			G(i,j)=GAUDEV(W_mean, W_sig, iseed)*G0(i,j)
			! G(j,i)=gau_mean*G0(j,i)
			G(j,i)=GAUDEV(W_mean, W_sig, iseed)*G0(j,i)
			enddo; enddo
		endif
	end if

    do i = 1, n
		indegree(i) = sum(G0(i,:))
		outdegree(i) = sum(G0(:,i))
        inweightdegree(i) = sum(G(i,:))
        outweightdegree(i) = sum(G(:,i))
	enddo	
    nedge1 = sum(G0)
    print*, "generate network done"
	write(6,*) repeat("=", 40)
    
    r = 0d0; noise = 0d0
    do i = 1, n
	   r(i) = r0*(ran3(iseed)+.5)
!	   a(i) = i*1d0
	enddo
    noise_diag = 0d0; noise_corr = 0d0
    do i = 1, n
	   noise_diag(i,i) = sqrt(DD)
       noise_corr(i,i) = 1d0
	enddo
	
    do i = 1, n-1; do j = i+1, n
       noise_corr(i,j) = 0.5d0
       noise_corr(j,i) = noise_corr(i,j)
    enddo; enddo

    noise = matmul(noise_diag, matmul(noise_corr, noise_diag))	
    noise_cholesky = noise
    call DPOTRF('L', N, noise_cholesky, N, info) ! Cholesky decomposition: noise = L L^T
    if (info .ne. 0) then
        print*, "Error in Cholesky decomposition, info = ", info
        stop
    endif
    ! 只保留下三角，避免上三角殘值干擾
    do i = 1, N
        do j = i+1, N
            noise_cholesky(i,j) = 0d0
        enddo
    enddo

    do i = 1, N;
        write(11,*) i, r(i)
    do j = 1, N
        write(12,*) i, j, noise(i,j)
        if (G(i,j) .eq. 0d0) cycle
        write(13,*) i, j, G(i,j)
    enddo; enddo
    
    avex = 0d0; avef = 0d0; Ko = 0d0; Ktau = 0d0; KtauBUF = 0d0
    i_avg=0; ITIME=0; Icount=1; jcount=0; IT=0; JTEM=0
    do i = 1, N; x1(i) = ran3(iseed)+.5d0; f1(i) = 0d0; enddo ! initial condition

    call cpu_time(time1)	
    
    do ii = 1, NSTEP
        tt = t0 + ii*dt
        if (mod(ii, int(nstep*1./10)) == 1) then
	   	    call progress(ii, nstep)
	   	endif        
        do jj = 1, N; gau(jj) = GAUDEV(0.d0,1.d0,ISEED); enddo ! generate standard Gaussian variates for noise
        call Euler_Maruyama(N, r, G, f1, x1, gau, noise_cholesky, noise_increment)
        ! write(10,*) ii
        if (mod(ii, nsample) .eq. 0) then
            i_avg = i_avg + 1
            avex = avex + x1
            avef = avef + f1
            DO ia=1,N; DO ib=1,N            	
                Ko(ia,ib)=Ko(ia,ib)+x1(ia)*x1(ib)
            ENDDO; ENDDO
            IF(ITIME .LT. MAXT_Ktau) THEN
                do ib=1,N;Ktaubuf(ib,ITIME)=X1(ib);enddo
	  	    ELSE
                DO IT=1,MAXT_Ktau
        	        JTEM=MOD(ITIME-IT,MAXT_Ktau)
	       	        do ia=1,N;do ib=1,N	
                    Ktau(ia,ib,IT)=Ktau(ia,ib,IT)+X1(ia)*Ktaubuf(ib,jtem)
	       	        enddo;enddo
	     	    ENDDO	
	     	    JTEM=MOD(ITIME,MAXT_Ktau)	
                do ib=1,N; Ktaubuf(ib,jtem)=X1(ib) ;enddo
	            jcount=jcount+1	     
	  	    ENDIF
        	ITIME=ITIME+1
        ELSE	
          	Icount=icount+1
        END IF !nsample
    enddo
    avex=avex/DBLE(i_avg); avef=avef/DBLE(i_avg)
    Ko=Ko/DBLE(i_avg); Ktau=Ktau/(DBLE(i_avg)-MAXT_Ktau)

    do ia=1,N; do ib=1,N
	   	Ko(ia,ib)=Ko(ia,ib)-avex(ia)*avex(ib) !covariance matrix
	   	do it=1,MAXT_Ktau
	      	Ktau(ia,ib,it)=Ktau(ia,ib,it)-avex(ia)*avex(ib) !Ktau matrix
	   	enddo
	enddo; enddo

    call cpu_time(time2)

    write(*,*) "Done !"	
	write(*,*) "time (steady_state) = ", time2-time1
	write(6,*) repeat("=", 40)

    do i = 1, N
        write(14,*) i, avex(i)
        write(15,*) i, avef(i)
        do j = 1, N            
            write(17,*) i, j, Ktau(i,j,1)         
            if (i > j) cycle
            write(16,*) i, j, Ko(i,j)            
        enddo
    enddo

    call build_jacobian_from_state(N, r, G, avex, Q)

    ! calculate the Schur form of Q
    call schur_factorize_real_matrix(N, Q, QschurT, QschurVS, QschurWR, QschurWI, QschurWORK, QschurBWORK)
    call solve_qx_xqt_from_schur_real(N, QschurT, QschurVS, -noise, QschurSol)
	KoT = cmplx(QschurSol, 0d0, kind=8)

    Lyapunov_check = matmul(Q, Ko) + matmul(Ko, transpose(Q)) + noise
    do i = 1, N; do j = 1, N
        write(18,*) Lyapunov_check(i,j)
    enddo; enddo
    do i = 1, N
		write(19,*) Ko(i,i), real(KoT(i,i), kind=8), aimag(KoT(i,i))
		do j = i+1, N		
        	write(20,*) Ko(i,j), real(KoT(i,j), kind=8), aimag(KoT(i,j))
    enddo; enddo

    end program Correlated_Noise

    subroutine Euler_Maruyama(n, r, G, f1, x1, gau, noise_cholesky, noise_increment)
	use global_vars
	implicit double precision(a-h, o-z)
	integer :: n
	real(8), intent(in) :: G(n,n), r(n), gau(n), noise_cholesky(n,n) 
    real(8) :: noise_increment(n)
	real(8) :: x1(n), f1(n), Gi(n)
	! One Euler-Maruyama update:
	!   x(t+dt) = x(t) + [r*x*(1-x) + Gi] dt + gau*sqrt(dt)
	!   Gi(i)   = sum_j G(i,j)*(x_j - x_i)	

    noise_increment = matmul(noise_cholesky, gau) * sqrt(dt) ! 產生 correlated noise 的增量

	Gi=0
    do jn=1,N    
        sumg=0.d0
        do j=1,N
            sumg=sumg+G(jn,j)*(X1(j)-X1(jn)) !(x2-x1)               
        enddo
        Gi(jn)=sumg
    enddo
!	f1 = r*(a-x1) + Gi
	f1 = r*x1*(1-x1) + Gi	
	x1 = x1 + f1 * dt + noise_increment ! x(t_{i+1})
	end subroutine Euler_Maruyama

    subroutine build_jacobian_from_state(N, r, G, xref, Q)
	implicit none
	integer, intent(in) :: N
	real(8), intent(in) :: r(N), G(N,N), xref(N)
	real(8), intent(out) :: Q(N,N)
	integer :: i
	real(8) :: inweightdegree_i
	! Logistic local dynamics:
	!   f_i(x_i) = r_i x_i (1 - x_i)
	! gives f_i'(x_i) = r_i (1 - 2 x_i).
	! Diffusive coupling h(x_i, x_j) = x_j - x_i contributes:
	!   off-diagonal:  +G(i,j)
	!   diagonal:      -sum_j G(i,j)
	Q = G
	do i = 1, N
		inweightdegree_i = sum(G(i,:))
		Q(i,i) = Q(i,i) + r(i) * (1d0 - 2d0*xref(i)) - inweightdegree_i
	enddo
	end subroutine build_jacobian_from_state

    logical function select_eig(wr, wi)
	implicit none
	real(8), intent(in) :: wr, wi
	! Callback for DGEES; no reordering needed in this code path.
	select_eig = .false.
	end function select_eig

    subroutine schur_factorize_real_matrix(N, Q, Tschur, VS, WR, WI, WORK, BWORK)
		implicit none
	integer, intent(in) :: N
	integer :: INFO, SDIM
	real(8), intent(in) :: Q(N,N)
	real(8), intent(out) :: Tschur(N,N), VS(N,N), WR(N), WI(N)
	real(8), intent(inout) :: WORK(*)
	logical, intent(inout) :: BWORK(N)
	logical, external :: select_eig

	Tschur = Q
	call DGEES('V', 'N', select_eig, N, Tschur, N, SDIM, WR, WI, VS, N, WORK, 8*N, BWORK, INFO)
	if (INFO .ne. 0) then
		write(6,*) "Warning: DGEES failed in schur_factorize_real_matrix, INFO =", INFO
	endif
	end subroutine schur_factorize_real_matrix

    subroutine solve_qx_xqt_from_schur_real(N, Tschur, VS, RHS, X)
		implicit none
		integer, intent(in) :: N
		real(8), intent(in) :: Tschur(N,N), VS(N,N), RHS(N,N)
		real(8), intent(out) :: X(N,N)
		real(8), allocatable :: TMPR(:,:), C(:,:)

		allocate(TMPR(N,N), C(N,N))
		call solve_qx_xqt_from_schur_real_ws(N, Tschur, VS, RHS, TMPR, C, X)
		deallocate(TMPR, C)
	end subroutine solve_qx_xqt_from_schur_real

    subroutine solve_qx_xqt_from_schur_real_ws(N, Tschur, VS, RHS, TMPR, C, X)
		! Solves the Sylvester equation Q*X + X*Q^T = RHS for X, where Q is given in real Schur form as Q = VS*Tschur*VS^T.
		implicit none
		integer, intent(in) :: N
		integer :: i, j, INFO
		real(8), intent(in) :: Tschur(N,N), VS(N,N), RHS(N,N)
		real(8), intent(inout) :: TMPR(N,N), C(N,N)
		real(8), intent(out) :: X(N,N)
		real(8) :: SCALE, symv
	
		! Solve Q*X + X*Q^T = RHS using the cached real Schur form Q = VS*Tschur*VS^T.
		call dgemm('N', 'N', N, N, N, 1d0, RHS, N, VS, N, 0d0, TMPR, N)
		call dgemm('T', 'N', N, N, N, 1d0, VS, N, TMPR, N, 0d0, C, N)
		call DTRSYL('N', 'T', 1, N, N, Tschur, N, Tschur, N, C, N, SCALE, INFO)
			if (INFO .ne. 0) then
				write(6,*) "Warning: DTRSYL returned INFO in solve_qx_xqt_from_schur_real_ws =", INFO
			endif
		if (SCALE .gt. 0d0 .and. abs(SCALE-1d0) .gt. 1d-14) C = C / SCALE
		call dgemm('N', 'T', N, N, N, 1d0, C, N, VS, N, 0d0, TMPR, N)
		call dgemm('N', 'N', N, N, N, 1d0, VS, N, TMPR, N, 0d0, X, N)
	
		do i = 1, N
			do j = i+1, N
				symv = 0.5d0*(X(i,j) + X(j,i))
				X(i,j) = symv
				X(j,i) = symv
			enddo
		enddo
	end subroutine solve_qx_xqt_from_schur_real_ws

    subroutine progress(i, nprog)
	implicit double precision(a-h,o-z)
	integer, intent(in) :: i, nprog
	! Simple single-line progress output.
		
	write(*,*) 
	write(*, '(a17)', advance='no') ' Progress : 10.0%'
	
	val = i*1./nprog+0.1
	
	write(*,'(a, a, f5.1, a2, a1)', advance="no") char(13), ' progress: ', val*100, '%'		
	
	end subroutine progress

    FUNCTION GAUDEV(RMU,SIG,ISEED)  ! Gaussian random numer generator
        implicit double precision(a-h,o-z)
	! Box-Muller transform with uniform variates from RAN3.
        TPI=6.28315307d0
	rand1=RAN3(iseed)
	DO WHILE (rand1==0d0)
	 rand1=RAN3(iseed)
	END DO

	rand2=RAN3(iseed)
	DO WHILE (rand2==0d0)
	 rand2=RAN3(iseed)
	END DO

       GAUDEV=RMU+SIG*SQRT(-2.d0*LOG(rand1))*COS(TPI*rand2)

        RETURN
    END


    FUNCTION RAN3(IDUM)
        INTEGER IDUM
        INTEGER MBIG,MSEED,MZ
        REAL(kind=8) RAN3,FAC
!       Subtractive RNG (Numerical Recipes style), output in [0,1).
!         PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=2.5E-7)
!   the routine returns random number in [0,1)  cccccccccccc
       PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=.999999E-9)

!     FAC=.999999E-9)
        DIMENSION MA(55)
        SAVE IFF,INEXT,INEXTP,MA
       DATA IFF /0/
       IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1

        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE

        DO 13 K=1,4
          DO 12 I=1,55

            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG

12        CONTINUE
13      CONTINUE

        INEXT=0
        INEXTP=31
        IDUM=1
      ENDIF

      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      RETURN
      END
