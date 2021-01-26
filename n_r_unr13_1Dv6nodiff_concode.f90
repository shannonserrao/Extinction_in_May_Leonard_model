!------------------------------------------------------------------------------------------------------------------------------------------
module all_in ! saving time evolution in seperate files

	!use random_gen ! module that contains Poisson distribution generator

	implicit none
	! global varibles:
	! N  			- number of species
	! p  			- number of predators (and prays - if not symmetric p=0)
	! nx 			- number of rows
	! ny 			- number of columns
	! n_iteration	- number of iterations unless all but one specie go extinct
	! ws 			- writing step
	! seed 			- seed for random number generator
	! mu_ 			- lower (_min) and upper(_max) boundaries to draw mean of PD
	! V				- volume of a lattice site = max(V,{\sum_i(s(i))})
	! latt_size		- nx*ny

	type time
		real*8	:: rt
		integer*4	:: ls
	end type time

	integer*4 				:: N, nx, ny, ws, seed, pray, mu_min, mu_max, V, latt_size, nx2, middle, my, mx, n_iteration, m, nn
	integer*4				:: fp_a, fp_b, fp_c, timecheck, tot_time
	!character				:: fm*17, fname*30, path*48, path2*63, csd*3
	character				:: fm*17, fname*10, path*35, path2*46, csd*3, jobfolder*14
	real*8					:: zero = 10.0**(-33), rate_d, rate_k, rate_p, rate_r
        double precision                        :: tester
	logical 				isext
	type(time)				ext_time_store
	contains
!_______________________________________________________________________

		subroutine calculations

		type(time), allocatable	:: tau(:)
		type(time)				:: help
		real*8, allocatable		:: k(:,:), r(:), d(:), p(:), total(:)
		real*8, allocatable		:: rateK(:,:,:), rateR(:,:), rateD(:,:), rateP(:,:)
		integer*4, allocatable	:: seed_array(:), s(:,:), ext_array(:), ext(:) !totpop change
		real*8 					:: mu, rn, rate, rs, t, tot_old,tempvar
		integer*4					:: i, j, ns, si, sj, un, iter, mm, place, totext,iter2,  timevar
		integer*4					:: sn,  n_site, nn_site,  ext_sum, s_tot, site, datetime(8)
		logical					:: exL, diff, again, ext1, newext
		logical, allocatable	:: empty(:), extflags(:)
		character				:: cs*1, ca*1, ctime*5, csi*3,csp*1
		integer*4, allocatable  :: 	 totpop(:)
		! ext  flags for extinction
 		integer*4, allocatable                  :: extinction(:)
		!

			! ext change set ext flag to zero and total extinctions to zero
			ext1 = .false.
			newext = .false.
            totext = 0
			! ext change
			allocate(k(1:N,1:N), r(1:N), d(1:N), p(1:N))

			k = 0.d0
			do i = 1, N
				do j = 1, pray
					site = mod(i+j-1,n) + 1
					k(i,site) = rate_k
				end do
				r(i) = rate_r
				d(i) = rate_d
				p(i) = rate_p
			end do

			  call date_and_time(VALUES=datetime)
                          seed = datetime(4) * (360000*datetime(5) + 6000*datetime(6) + 100*datetime(7) + datetime(8)) + seed
				seed = abs(seed)
				!print *, seed

				!print *, datetime(4)
			!version 4 seed program
			! set seed for the random number generator
			call random_seed(size = sn)
			allocate(seed_array(sn))
			seed_array = seed
			call random_seed(put=seed_array)
			deallocate(seed_array)

			! ext allocate extinction flag
			allocate(s(1:N,1:latt_size), tau(1:latt_size), total(1:latt_size), empty(1:latt_size))
		 	allocate(extinction(1:N), extflags(1:N), totpop(1:N)) !totpop(1:N)


			! set initial conditions - Poisson distribution
			! draw a random number from interval [µ_min, µ_max] for
			! each species and each site and then call Poisson distr. generator
			! with a chosen µ_{\alpha,i}, where \alpha is a species, i is a site

			do j = 1, latt_size
				do i = 1, N
					call random_number(rn)
					tempvar=rate_r/(rate_p+rate_k)
					if (isext .eqv. .false.) then
						s(i,j)= int(V*(tempvar)+1)
					endif
					totpop(i)=0
				end do
			end do



 			do i= 1, N
   				extinction(i) = 0
				extflags(i) = .false.
			end do





			k = k/dfloat(V)
			p = p/dfloat(V)

			allocate(rateR(1:N,1:latt_size), rateK(1:N,1:N,1:latt_size), rateD(1:N,1:latt_size),rateP(1:N,1:latt_size))

			if  (isext .eqv. .true.) then
				s(1,1)= fp_a
				s(2,1)=fp_b
				s(3,1)=fp_c
				!tau(1)%ls = 1
				!tau(1)%rt = ext_time_store
				!goto 100
			endif



			! initialaze times for the first (re)action at each site
			! rates are proportional to the number of species at the site
			do i = 1, latt_size
				call rates(i, k, r, d, p, s, total(i), rateR(1:N,i), rateK(1:N,1:N,i), rateD(1:N,i), rateP(1:N,i))
				call random_number(rn)
				rn = rn + zero
                    		tau(i)%ls = i
				tau(i)%rt = -log(rn)/total(i)
				if (isext .eqv. .true.) then
					tau(i)%rt = ext_time_store%rt - log(rn)/total(i)
				endif
			        !   if added
				  !if (isnan(tau(i)%rt)) then
					!		write (3,*) " Error in initial"
                                        !                write (3,*) log(rn)
                                        !                write (3,*) total(i)
					!		 call exit(0)
				  !end if


				!if (total(i) == 0 .or. rn == 0) then
				!	write(1,'(a)') 'Infinitity or Nan occured'
				!end if
			end do


			!100 !line 100

			! sort times from smallest to largest

			do i = 1, latt_size-1
				do j = i, latt_size
					if (tau(i)%rt > tau(j)%rt) then
						help = tau(j)
						tau(j) = tau(i)
						tau(i) = help
					end if
				end do
			end do




			allocate(ext_array(1:N), ext(1:N))

			!empty = .false.
			!ext_array = 0

			do mm = 1, 2**m

				! open files to store data, one file per site and per species

				do i = 1, latt_size
					do j = 1, N
						write(csp,'(i1.1)') j
						write(csi,'(i3.3)') i
						ctime = csp//'_'//csi
						write(csd,'(i3.3)') mm

					end do
				end do


					if (mm == 1) then
						do i = 1, latt_size
							do j = 1, N

								un = (i-1)*N+j+4
								totpop(j) = s(j,i) + totpop(j)
							end do
						end do
					! ext changes ----------------------------------
						do j=1,N
						  extinction(j) = totpop(j)
						  if (extinction(j)  == 0) then
							extflags(j) = .true.
							ext1 = .true.
							totext = totext + 1
				                  end if
					!ext changes	extinction(j) = totpop(j)

						  totpop(j)=0
						end do
						!write(3,*) tau(1)%rt

   					          !ext changes
						if ((ext1 .eqv. .true.) .and. (totext == N-1)) then
						  !write(4, '(i11.10)  ', advance ='no') seed
						  write(4, '(i10.10) ', advance ='no') seed
						  write(4, '(a) ', advance ='no') '  '
						  write(4,'  (i2.2) ', advance ='no') totext
						  do j=1,N
 							write(4,'(i6.5) ', advance ='no') extinction(j)
						  end do
							write(4,*) tau(1)%rt
						  ext1 = .false.
                                                end if

					! ext changes





						!write(2,*) tau(1)%rt
					end if

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! start iteration
				do iter2 = 1, 1000
				do iter = 1, n_iteration 	! start simulation
				         ! write (3,*) tau(1)%rt
					t  = tau(1)%rt

                                        site = tau(1)%ls		! (re)act on a site with the smallest waiting time tau(1)
					exL  = .false. 			! will be set to true when at least one of the species at some time is extinct
					rate = zero





						! choose random number to find a reaction
						call random_number(rn)


						do i = 1, N
							! if not diffusion check if anihilation
							rate = rate + rateP(i,site)
							if ( (rate - rn) > zero ) then
								s(i,site) = s(i,site) - 1
								ca = 'p'
								goto 200
							end if
						end do


						! if not anihilation check if predation
						do i = 1, N
							do j = 1, N
								rate = rate + rateK(j,i,site)
								if ( (rate-rn) > zero ) then
									s(i,site) = s(i,site) - 1
									ca = 'k'
									goto 200
								end if
							end do
						end do

						do i = 1, N
							! if not predation check if reproduction
							rate = rate + rateR(i,site)
							if ( (rate - rn) > zero ) then
								s(i,site) = s(i,site) + 1
								ca = 'r'
								goto 200
							end if
						end do !timee(8)



					goto 100 				! nothing happened - calculate new tau(site)



					! no diffusion so 300 should not be exec ..
					! if diffusion recalculate the rates on the n_site site (number of individuals changed)
					300    tot_old = total(n_site)
                    				call rates(n_site, k, r, d, p, s, total(n_site), rateR(1:N,n_site), &
                                    		&  rateK(1:N,1:N,n_site), rateD(1:N,n_site), rateP(1:N,n_site))
                    				do i = 1, latt_size
                        				if (tau(i)%ls == n_site) then
                          				  si = i
                         				   exit
                       				         end if
                    				end do


                 				   400 again = .false.
                 				   call random_number(rn)
                  	   		           rn = rn + zero
                                                   help%rt = -log(rn)/total(n_site) + t
                                                   help%ls = n_site

                    				   call queue(help, si, tau, again)
                                                   if (again) goto 400



					! recalculate rates at site
					200 call rates(site, k, r, d, p, s, total(site), rateR(1:N,site), rateK(1:N,1:N,site), rateD(1:N,site), rateP(1:N,site))

					! for no diffusion or nothing happened calculate new tau just at site, and put it in the queue
					100  call random_number(rn)
						rn = rn + zero

						if (total(site) < zero) then
							help%rt = tau(latt_size)%rt *(1.0 + real(1/latt_size))
							help%ls = site
							again= .true.
						else
							! write (3,*) 'im in 100'
							help%rt = -log(rn)/total(site) + t
							help%ls = site
							again = .false.
							!si = 1
						end if
						!call queue(help, 1, tau, again) ! changes made for diffusion

						tau(1)%rt = help%rt
						tau(1)%ls = help%ls
						!end changes made
						if (again) goto 100


						! write number of species on sites
						do i = 1, latt_size
							do j = 1, N
								! ext change j+3 to j+4
								!un = (i-1)*N+j+4
								!write(un) s(j,i)
								!write(un,'(i10)') s(j,i)
								totpop(j)= s(j,i) + totpop(j)
							end do

						end do

                                        !ext changes
						do j=1,N
						   extinction(j) = totpop(j)
						  if (extinction(j)  == 0) then
							if (extflags(j) .eqv. .false.) then
								totext = totext + 1
								newext = .true.
							end if
							extflags(j) = .true.
							!ext1 == .true.
				                  end if
					!ext changes 	totpop(j) = extinction(j)
						  !write(3,'(i6.5) ', advance ='no') totpop(j)
						  totpop(j)=0
						end do
						!write(3,*) tau(1)%rt
				       !test-----
						     !ext changes
						if ((newext .eqv. .true.)) then ! .and. (totext == N-1)) then
						  write(4, '(i10.10)  ', advance ='no') seed
						  write(4, '(a) ', advance ='no') '  '
						  write(4,'(i2.2) ', advance ='no') totext
						  do j=1,N
 							write(4,'(i6.5) ', advance ='no') extinction(j)
						  end do
							write(4,*) tau(1)%rt
						  newext = .false.
                                                end if
						  ! ext changes

						if (totext == (N-1)) exit


					!end if
					! concode
					call system_clock(timevar)
					if ((timecheck - timevar) < 1000) then
						if ((totpop(1)+totpop(2)+totpop(3)) > 0) then
							open(unit=6,file=path2//'_'//'ext.dat',form = 'formatted', access='sequential')
							do j=1,N
								write(6,'(i6.5) ', advance ='no' ) totpop(j)
							end do
							write(6,*) tau(1)%rt
							close(6)
							print *, 'WARNING: Extinction not reached!'
						end if
					end if

					!concode
				end do ! iteration
					if (totext == (N-1)) exit
				end do ! iter2
				print *, 'last time recorded is-', tau(1)%rt
				if (totext /= (N-1)) then
						print *, 'The number of extinctions is -', totext
						print *, 'the total populations are - '
						do j = 1, N
				    			print *, totpop(j)
						end do ! pop print
				end if
				exit
				! close all open files (1=parameters, 2=time, else=data)
				! ext change +2 to +3
				!do i = 2, N*latt_size+
				!	close(i)
				!end do
				!close(4)


			end do !mm

			! at the end of simulation write ext_array to parameters file
			!write(1,'(a)') 'Extinction times, if 0, no extinction'
			!write(1,*) ext_array


			!close(3)
			deallocate(k, r, d, p, s, ext_array, ext, total, tau, empty,totpop,extinction,extflags)
			deallocate(rateR, rateK, rateD, rateP)


		end subroutine
!_______________________________________________________________________

		subroutine queue(help, takeout, tau, again)

			type(time), intent(inout)	:: tau(1:latt_size)
			type(time)					:: tau_n(1:latt_size)
			type(time), intent(in)		:: help
			logical, intent(inout)		:: again
			integer*4, intent(in)			:: takeout
			integer*4						:: i, j, g0, g1, g2, h0, h1 , place


			if (takeout == 1) then
				do i = 1, latt_size-1
					tau_n(i) = tau(i+1)
				end do
			elseif (takeout == latt_size) then
				do i = 1, latt_size-1
					tau_n(i) = tau(i)
				end do
			else
				if (takeout > latt_size) then
				   Print *, 'error in the takeout number'
				end if
				do i = 1, takeout-1
					tau_n(i) = tau(i)
				end do
				do i = takeout, latt_size-1
					tau_n(i) = tau(i+1)
				end do
			end if

			tau_n(latt_size)%rt = tau(latt_size)%rt + help%rt
			tau_n(latt_size)%ls = 0

			! check if it is smaller than the first one

			if ((help%rt - tau_n(1)%rt) < zero) then
				place = 1
				goto 2000
			end if



			g1 	= latt_size / 2
			h1 	= g1 + 1
			g0  = g1


			do i = 1, nn - 1

				if (abs(help%rt - tau_n(h1)%rt) < zero) then
					! same tau, recalculate help
					again = .true.
					goto 1000

				elseif (help%rt - tau_n(h1)%rt < zero) then

					g1  = g0 - 2**(nn-1-i)
					h1  = g1 + 1
					g0  = g1
				        if ( (h1-1) > latt_size) then
						print *, 'g1- ', g1, 'h1- ', h1, 'g0- ', g0
					end if
				else

					g1  = g0 + 2**(nn-1-i)
					h1  = g1 + 1
					g0 = g1
					!print *, 'g1- ', g1, 'h1- ', h1, 'g0- ', g0
					if ( (h1-1) > latt_size) then
						print *, 'g1- ', g1, 'h1- ', h1, 'g0- ', g0
					end if

				end if
			end do
			        if (h1 > latt_size) then
				 Print *, 'there is an error in h1'
				end if
				if (help%rt - tau_n(h1)%rt < zero) then
					place = h1
				else
					place = h1 + 1
				end if

				if (place == latt_size) goto 3000

				2000 do i = latt_size , place+1, -1
					tau_n(i) = tau_n(i-1)
				end do


				3000 tau_n(place) = help

				tau = tau_n

				1000 return

		end subroutine

!_______________________________________________________________________

		subroutine neighbor(site,dir,n_site)

			integer*4, intent(in)	:: site, dir
			integer*4				:: si, sj, n_site


			select case (dir)
				case (1) ! right neighbor
					n_site = mod(site, ny) + 1
					!print *, 'right ', site, n_site
				case (2) ! left neighbor
					n_site = ny - mod (ny - site +1, ny)
					!print *, 'left ', site, n_site
			end select

		end subroutine

!_______________________________________________________________________

		subroutine rates(site, k, r, d, p, s, tot, rR, rK, rD, rP)

			real*8, intent(in)	:: k(1:N,1:N), r(1:N), d(1:N), p(1:N)
			integer*4, intent(in)	:: site, s(1:N,1:latt_size)
			real*8, intent(out)	:: rK(1:N,1:N), rR(1:N), rD(1:N), rP(1:N)
			real*8, intent(out)	:: tot
			integer*4				:: i, j, n_site, ns

			tot = zero
			do i = 1, N
				do j = 1, N
					tot = tot + (k(j,i)*s(j,site))*s(i,site)
				end do
				tot = tot + r(i)*s(i,site) +  2.d0*d(i)*s(i,site) + p(i)*s(i,site)*(s(i,site)-1)

			end do

			if (tot < zero) then
					rK = zero
					rR = zero
					rP = zero
					rD = zero

			else
				do i = 1, N
					do j = 1, N
						rK(j,i) = k(j,i)*s(j,site)*s(i,site)/tot
					end do

					rR(i) = r(i)*s(i,site)/tot
					rP(i) = p(i)*s(i,site)*(s(i,site)-1)/tot
					rD(i) = d(i)*s(i,site)/tot

				end do
			end if





		end subroutine

!_______________________________________________________________________



end module

!_______________________________________________________________________


! main program

program NspeciesGame
use all_in
!include 'mpif.h'
implicit none

include 'mpif.h'
         ! version 5 change statistics iterations
          integer*4       :: ind, maxind, qsubid
         !
	character	:: arg*30, cn*1, cp*1, cx*4, cy*4, da*8, ti*10, broj1*8, broj2*6,cqid*7,cpid*2,procfile*8,Vol*6,kap*4,pfname*12
	integer*4		:: beginning, r_ate, can, en, i, j, timee(8)
	integer*4         ::  ierr, num_procs, procid

	!logical :: file_exists
	call system_clock (beginning, r_ate)
	call DATE_AND_TIME(values=timee)
	call SYSTEM_CLOCK(timecheck)
        print *, timee
	print *, 'The time now is ', timecheck
	timecheck=timecheck+tot_time
	PRINT *, 'Time before the job abends is ', timecheck
        maxind=200
	num_procs=1
	isext=.false.
	! default values: (can be changed from the command line)

	do j = 1, 1

		N 				= 3
		pray			= 1 ! if we choose non-symmetrical game p = 0
		nx 				= 1
		nn				= 0
		ny 				= 2**nn
		n_iteration 	= 2**30
		m				= 0
		ws 				= 2**15
		seed			= 1
		mu_min			= 15
		mu_max			= 25
		V				= 10
		rate_p			= 1.0  ! p is gamma the death rate
		rate_r			= 0.1   ! rho or birth rate
		!rate_p			= 10.d0
		!rate_r			= 1.d0
		qsubid=1
		procid=1

		select case (j)
			case (1)
				rate_d			= 0.0  !/dfloat(ny**2)
				rate_k			= 1.00  ! the predation rate
			case(2)
				rate_d			= 1.d0   !/dfloat(ny**2)
				rate_k			= 1.0d0
			case(3)
				rate_d			= 10.d0   !/dfloat(ny**2)
				rate_k			= 1.0d0
			case (4)
				rate_d			= 0.1d0   !/dfloat(ny**2)
				rate_k			= 1.0d0
			case(5)
				rate_d			= 1.d0   !/dfloat(ny**2)
				rate_k			= 1.0d0
			case(6)
				rate_d			= 10.d0   !/dfloat(ny**2)
				rate_k			= 1.0d0
		end select

		can = command_argument_count()


                do i = 1, can-1, 2
			call getarg(i,arg)
				if 	   ( arg == "nx" ) then
					call getarg(i+1,arg)
					write(nx,'(i3)')arg
				elseif ( arg == "nn" ) then
					call getarg(i+1,arg)
					read(arg,'(i3)') nn
				elseif (arg == "ni") then
					call getarg(i+1,arg)
					write(n_iteration,'(i12)') arg
				elseif (arg == "ws") then
					call getarg(i+1,arg)
					write(ws,'(i12)') arg
				elseif (arg == "N") then
					call getarg(i+1,arg)
					write(N,'(i1)') arg
				elseif (arg == "p") then
					call getarg(i+1,arg)
					write(pray,'(i1)') arg
				elseif (arg == "si") then
					call getarg(i+1,arg)
					write(seed,'(i5)') arg
				elseif (arg == "qd") then
					call getarg(i+1,arg)
					read(arg, '(i7)') qsubid
				elseif (arg== "pd") then
					call getarg(i+1,arg)
					read(arg, '(i2)') procid
				elseif (arg== "V") then
					call getarg(i+1,arg)
					read(arg, '(i6)') V
					!read(arg, *) Vol
				elseif (arg== "r") then
					call getarg(i+1,arg)
					read(arg,*) rate_r
				elseif (arg== "g") then
					call getarg(i+1,arg)
 					read(arg,*) rate_p
				elseif (arg== "k") then
					call getarg(i+1,arg)
 					read(arg,*) rate_k
					read(arg,*) kap
				elseif (arg== "d") then
					call getarg(i+1,arg)
					read(arg,*) rate_d
				elseif (arg== "mi") then
					call getarg(i+1,arg)
					read(arg, '(i7)') maxind
				elseif (arg== "EP") then
					call getarg(i+1,arg)
					read(arg, *) jobfolder
					isext=.true.
				elseif (arg== "TT") then
					call getarg(i+1,arg)
					read(arg, '(i8)') tot_time

				endif
		enddo
		! file changes here
		!goto 5000
		!Mpi begins
                call MPI_INIT (ierr)

		call MPI_COMM_RANK (MPI_COMM_WORLD, procid, ierr)
		write(cpid,'(i2.2)') procid
		procfile= cpid//'_input.dat'

	        ny=2**nn
		write(Vol, '(i6.6)')V
		write(cx,'(i4.4)') nx
		write(cy,'(i4.4)') ny
		write(cn,'(i1.1)') N
		write(cp,'(i1.1)') pray
		call date_and_time(da,ti)
		if (isext .eqv. .false.) then
			read(da,'(A8)'),broj1
			read(ti,'(A6)'),broj2
		end if
		write(cqid,'(i7.7)') qsubid
		write(cpid,'(i2.2)') procid

		!path  = '/media/dl/external-DL/Research/GameTh/Num/Data/'
		!path  = '/home/dl/Research/GameTh/Num/Data/'//cn//'_'//cp//'/'//cx//'x'//cy//'/'
		path  = 'shann87/fortran/Data/'//cn//'_'//cp//'/'//cx//'x'//cy//'/'
                !f	name = '('//cn//'_'//cp//')_'//cx//'x'//cy//'_'//broj1//broj2




		do ind =1,5
			if (kap(ind:ind) == '.') then
				kap(ind:ind) = 'p'
			end if
		end do
		pfname = Vol//'_'//kap
		!print *,path2
		print *,pfname
		!open(unit=5,file=path//pfname//'.dat',form = 'formatted', access='sequential')
		!write(5, *) path2
		!close(5)

		!fname = qsubid//'_'//procid
		if (isext .eqv. .false.) then
			path2 = path//pfname//'/'
			print *, path2
		else
			path2 = path//jobfolder//'/'
		end if
		!print*,path
		!*  CHANGES
		!write(kap,'(f3.2)') rate_k
		!*	CHANGES
		!INQUIRE(FILE=path2//'_'//'ext.dat', EXIST=file_exists)----> changes ---> //broj1//broj2


		if (isext .eqv. .false.) then
			call system('mkdir -p '//path//pfname//'/' )
		else
			call system('mkdir -p '//path//jobfolder )
			open(unit=6,file=path2//'_'//'ext.dat',form = 'formatted', access='sequential')
			read (6 , * ) fp_a, fp_b, fp_c, ext_time_store
			close(6)
			!read (6 ,'(i6.5) ') fp_b
			!read (6 ,'(i6.5) ') fp_c
			!read (6, '(i6.5) ') ext_time_store
		end if


		! make m subfolders
		!do i = 1, 2**m
			!write(csd,'(i3.3)') i
			!print*,path
		        !print*,path2
			!call system('mkdir -p '//path2//csd )
			!call system('mkdir -p '//path2 )
		!end do
		print *, 'lattice dimension -' , nx, 'X', ny
		print *, 'volume -', V, ' gamma -' ,rate_p, ' rho -', rate_r, 'kappa -', rate_k, 'diff rate -', rate_d
		print *, 'N-R is ', N, '-', pray
		!print *, path
		!print *, path2


		!write(broj2,'(i6.6)') nx*ny
		!fm = '('//broj2//'(i1,2x))'

		latt_size = nx*ny
	        !open(unit=4,file=path2//csd//'/'//fname//'_'//'ext.dat',form = 'formatted', access='sequential')
		! MPI begins
		!call MPI_INIT (ierr)

		!call MPI_COMM_RANK (MPI_COMM_WORLD, procid, ierr)
		call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
		print *, "Hello world! I'm process ", procid, " out of ", num_procs, " processes."
		! MPI code ends
		write(cqid,'(i7.7)') qsubid
		write(cpid,'(i2.2)') procid
		fname = cqid//'_'//cpid
		!code
		!open(unit=5,file=path//'/'//pfname//'_'//'.dat',form = 'formatted', access='sequential')
		!write(5, *) path2
		!close(5)
		!*	CHANGES
		!
		open(unit=4,file=path2//'/'//fname//'_'//'ext.dat',form = 'formatted', access='sequential')
		do  ind= 1, maxind
		      seed= ind + 13*qsubid+ 17*procid
                      call calculations

	        end do
                close(4)
	end do
	! close extinction record file

	call DATE_AND_TIME(values=timee)
        print *, timee
	call system_clock(en)
	print*,real(en - beginning) / real(r_ate)
! MPI code begins
call MPI_FINALIZE (ierr)
!stop
!end
!5000 print *, 'hello'
!MPI code ends
end program
